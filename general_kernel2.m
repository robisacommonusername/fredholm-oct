%A better implementation of general_kernel, with a custom numerical
%integration routine. Better able to handle the 1/sqrt(k^2-q^2) term,
%due to incorporation of Gauss-Chebyshev quadrature

function H = general_kernel2(g0,Q,z0)
	s = @(QQ,rho,t,k) QQ-k*rho*[cos(t),sin(t)]; %Helper function, Q-q
	H = @(k,z) integrate_kernel_custom(...
		@(rho,t) k*rho*exp(i*k*(z-z0)*sqrt(1-rho^2))*g0(k*rho*[cos(t),sin(t)],k)*exp(i*k*(z-z0)*sqrt(1-norm(s(Q,rho,t,k))^2))*g0(s(Q,rho,t,k),k),1,0,2*pi,0.0001,10);
end

%Below is a custom adaptive quadrature method, for integrating a function
%with a radial weighting function over a disk or sector. The integrals 
%have the form
% \int_0^a \int_{t_a}^{t_b} \frac{f(r,t)}{\sqrt{a^2-r^2}} dt dr
%The radial integral uses a Gauss-Chebyshev quadrature method, and the
% angular part is by Gauss-Legendre quadrature.
%Don't include the radial weighting function in f, it will be included
%implicitly. For integral over range r \in [0,a), the implicit weighting
%function is 1/sqrt(a^2-r^2)
%Adaptive method: double number of radial points, and divide angular
%sub interval into two.
%For more details, see the implementation notes
function [I, status] = integrate_kernel_custom(f,r,ta,tb,varargin)
	tol = 0.0001; %absolute tolerance
	max_iters = 15; %maximum recursion depth =>max sub-sectors = 2^15 ~ 32k
	n=20; %initial num quad points for gauss-chebyshev
	if nargin > 4
		if varargin{1} > 0
			tol = varargin{1};
		end;
	end;
	if nargin > 5
		%Octave interpreter stack limits recursion depth to 255
		if varargin{2} < 255
			max_iters = varargin{2};
		else
			max_iters = 254;
		end;
	end;
	if nargin > 6
		n = varargin{3};
	end;
	if nargin > 7
		%previous estimate passed forward to prevent recalculation
		I1 = varargin{4};
	else
		I1 = 0;
	end;
	
	I=0;
	if max_iters > 0
		%split sector, and double number of radial points
		mt = (ta+tb)/2;
		I2a = integrate_chebyshev(f,r,ta,mt,2*n);
		I2b =integrate_chebyshev(f,r,mt,tb,2*n);
		I2 = I2a + I2b;

		if abs(I1-I2) < tol || abs(I2) < tol
			I = I2;
			status = 0;
		else
			%pass previous estimates forward into child calls to prevent
			%recalculation (i.e. pass I2a and I2b)
			[I1,stat1] = integrate_kernel_custom(f,r,ta,mt,tol/2,max_iters-1,2*n,I2a);
			[I2,stat2] = integrate_kernel_custom(f,r,mt,tb,tol/2,max_iters-1,2*n,I2b);
			status = max([stat1,stat2]);
			I = I1+I2;
		end;
	else
		%recursion limit reached
		I = I1;
		status = 1;
	end;
end;

%Non adaptive Gauss-Chebyshev/Gauss-Legendre integrator.
%See integrate_kernel_custom above and implementation notes.
%Also: http://dlmf.nist.gov/3.5#Ex7
function I = integrate_chebyshev(f,r,ta,tb,n)
	[pts_legendre,weights_legendre] = generate_quadrature('gauss10', 10);
	pts_t = pts_legendre*(tb-ta) + ta;
	weights_t = weights_legendre*(tb-ta);
	acosXkPlus = pi*(1:2:(n-1))/(2*n);
	pts_chebyshev = cos(acosXkPlus);
	weights_chebyshev = pi/n*ones(length(acosXkPlus),1);
	pts_r = r*pts_chebyshev;
	weights_r = weights_chebyshev;
	
	[rho,theta] = meshgrid(pts_r,pts_t);
	I = weights_t'*arrayfun(f,rho,theta)*weights_r;
end;
