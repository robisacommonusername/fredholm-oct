%This function is a prototype for implementation in C

%This function supersedes
%warp_variables, gauss_kernel, discretise_operator

%TODO:
%optimisations for when q<<k (linear phase approximation)
%Optimisation for rotational symmetry (function of q only => only need single integral)

function f = fastcall_general_kernel_disc(Qx, Qy, g0, q0,z0,varargin)
	if nargin > 5
		symmetric = varargin{1};
	else
		symmetric = 0;
	end;
	f = @(A,ki,zf,opts) fastcall_general_disc_worker(A,ki,zf,Qx,Qy,g0,q0,z0,symmetric,opts);
end

%This is the fastcall function, and it should be implemented in C
%One difficulty here is that g0 must be a function handle, which 
%introduces an enormous amount of interpreter overhead.

%Possibly we should just use the g0 and q0 as arguments to a code generator,
%which can then generate a mex/oct file which takes Qx and Qy as arguments
%however, for now we'll just write a pure matlab/octave implementation
function [Kd,Kdag,pts,weights] = fastcall_general_disc_worker(A,ki,zf,Qx,Qy,g0,q0,z0,symmetric,opts)
	if opts.low
		quad_method = opts.quad_method_low;
		n = opts.n_low;
	else
		quad_method = opts.quad_method;
		n = opts.n;
	end;
	
	ka = ki(1);
	kb = ki(end);
	if ka > kb
		ki = flipdim(ki);
		A = flipdim(A);
		ka = ki(1);
		kb = ki(end);
	end;
	
	[pts, weights] = generate_quadrature(quad_method, n);
	
	%Test for the non-overlapping case
	Q = sqrt(Qx^2 + Qy^2);
	if (Q >= 2*q0)
		Kd = zeros(n,n);
		Kdag = zeros(n,n);
	else
		if (abs(Q) < 10*eps)
			symmetric = 1;
			%Full overlap. This case needs to be handled specially
			%No change of variables needed
			height = 2*q0;
			width = 2*q0;
			centre = [0;0];
			cosPhi = 1;
			sinPhi = 0;
		else
			%partially Overlapping case - numerical integration
			
			%Find points where the boundary circles |qq| = q0 and |qq-QQ| = q0 intersect
			%^doubled indicates vector e.g. qq = (qx,qy)
			%Denote these intersection points as qc1 and qc2
			
			Qhat = 1/Q*[Qx;Qy];
			%The intersection points are found by rotating Qhat by an angle +- theta
			%and scaling the length to be q0
			cosTheta = Q/2/q0;
			sinTheta = sqrt(1-Q^2/4/q0^2);
			R = [cosTheta -1*sinTheta; sinTheta cosTheta];
			qc1 = q0*R*Qhat;
			qc2 = q0*R'*Qhat;
			
			%Now find bounding box for intersection region. Axes for box
			%are along line between qc1 and qc2, and perpendicular
			%bisector of this line segment.
			centre = 1/2*(qc1+qc2);
			height = norm(qc1-qc2);
			width = 2*q0-Q;
			phi = atan2(Qy,Qx); %angle between second box axis and new x axis
			cosPhi = cos(phi);
			sinPhi = sin(phi);
		end;
		%Coordinate transformation: translate centre to (0,0), then
		%rotate by -1*phi to transform q coordinates to new 
		%system of coords (x,y). Integrate over x in [-width/2,width/2]
		%y in [-height/2,height/2]. No jacobian factor required for
		%this transformation
		%Explicitly
		%[x;y] = [cosPhi -1*sinPhi;sinPhi cosPhi]*(q - centre)
		%[qx;qy] = [cosPhi sinPhi; -1*sinPhi cosPhi]*[x;y] + centre
		qx = @(x,y) cosPhi*x+sinPhi*y+centre(1);
		qy = @(x,y) cosPhi*y-sinPhi*x+centre(2);
		xa = -1/2*width;
		xb = 1/2*width;
		ya = -1/2*height;
		yb = 1/2*height;
		
		%resample A at the quadrature points. We don't use resample vector,
		%since we're going to need the k_quad points anyway, so there's no
		%point calculating them twice. Do it manually
		k_quad = unwarp_data('linear', pts, ka, kb);
		Abar = interp1(ki,A,k_quad,'linear','extrap');
		
		z_quad = unwarp_data(opts.warp_method, pts, zf);
		deriv = warp_derivative(opts.warp_method, pts, zf);
		
		%Perform the integration - test for some special cases that allow
		%us to optimise things
		status = 0;
		if q0/ka < 0.01
			%small q => linear phase approximation
			status = status + 1;
		end;
		if symmetric
			status = status+2;
		end;
		switch (status)
			case 1:
			%If q << k, then we can make the approximation sqrt(k^2-q^2) ~ k
			integrals = zeros(n,1);
			for ii = 1:n
				k = k_quad(ii);
				integrals(ii) = dbl_integrate_adaptive(@(x,y) 1/k*g0(qx(x,y),qy(x,y),k)*g0(Qx-qx(x,y),Qy-qy(x,y),k),...
					,xa,xb,ya,yb,0,7);
			end;
			[zz,kk] = meshgrid(z_quad,k_quad);
			K_unweighted = diag(integrals)*exp(2*i*kk*(zz-z0));
			
			case 2:
			%Symmetric case - no variable remapping is necessary, since region
			%is already a rectangle. This offers a substantial performance
			%boost, since annonymous functions are very slow.
			mid = norm(centre)/2;
			xa = mid-width/2;
			xb = mid+width/2;
			ya = -1*height/2;
			yb = 1*height;
			K_unweighted = zeros(n,n);
			ii = 1;
			for zi = 1:n
				z = z_quad(zi);
				for ki = 1:n
					k = k_quad(ki);
					%In the C version, the code generator should create a function
					%for the integrand expression which accepts centre, sinPhi, cosPhi, etc
					%such that the common sub-expressions only need to be computed once
					K_unweighted(ii) = dbl_integrate(...
						@(x,y) 2*pi*i/sqrt(k^2-x^2-y^2)*exp(i*(z-z0)*sqrt(k^2-x^2-y^2))*g0(x,y,k)*exp(i*(z-z0)*sqrt(k^2-(Qx-x)^2-(Qy-y)^2))*g0(Qx-x,Qy-y,k)...
						,xa,xb,ya,yb,0,7); %max recursion depth = 7 => no more than 1280 quad points/dimension
					ii = ii+1;
				end;
			end;
			case 3:
			%In this case, we use both optimisations, i.e. linear phase
			%And no variable remapping
			mid = norm(centre)/2;
			xa = mid-width/2;
			xb = mid+width/2;
			ya = -1*height/2;
			yb = 1*height;
			integrals = zeros(n,1);
			for ii = 1:n
				k = k_quad(ii);
				integrals(ii) = dbl_integrate_adaptive(@(x,y) 1/k*g0(x,y,k)*g0(Qx-x,Qy-y,k),...
					,xa,xb,ya,yb,0,7);
			end;
			[zz,kk] = meshgrid(z_quad,k_quad);
			K_unweighted = diag(integrals)*exp(2*i*kk*(zz-z0));
			otherwise:
			%Default integrator - slow
			%z increases across row, k decreases down column. Access in row
			%major order for efficiency. Bit of a moot point though, this
			%is going to be very slow in matlab/octave
			K_unweighted = zeros(n,n);
			ii = 1;
			for zi = 1:n
				z = z_quad(zi);
				for ki = 1:n
					k = k_quad(ki);
					%In the C version, the code generator should create a function
					%for the integrand expression which accepts centre, sinPhi, cosPhi, etc
					%such that the common sub-expressions only need to be computed once
					K_unweighted(ii) = dbl_integrate(...
						@(x,y) 2*pi*i/sqrt(k^2-(qx(x,y))^2-(qy(x,y))^2)*exp(i*(z-z0)*sqrt(k^2-(qx(x,y))^2-(qy(x,y))^2))*g0(qx(x,y),qy(x,y),k)*exp(i*(z-z0)*sqrt(k^2-(Qx-qx(x,y))^2-(Qy-qy(x,y))^2))*g0(Qx-qx(x,y),Qy-qy(x,y),k)...
						,xa,xb,ya,yb,0,7); %max recursion depth = 7 => no more than 1280 quad points/dimension
					ii = ii+1;
				end;
			end;
		end;
		
		K_unweighted = diag(Abar)*K_unweighted*diag(deriv);
		W = diag(weights);
		Kd = K_unweighted*W;
		Kdag = ctranspose(K_unweighted)*W;
	end;
	
end

