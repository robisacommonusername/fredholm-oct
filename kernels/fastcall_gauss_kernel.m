%This function is a prototype for implementation in C

%This function supersedes
%warp_variables, gauss_kernel, discretise_operator

function f = fastcall_gauss_kernel(Q,alpha,z0)
	f = @(A,ki,zf,opts) fastcall_gauss_worker(A,ki,zf,Q,alpha,z0,opts);
end

%This is the fastcall function, and it should be implemented in C
function [Kd,Kdag,pts,weights] = fastcall_gauss_worker(A,ki,zf,Q,alpha,z0,opts)
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
	%resample A at the quadrature points. We don't use resample vector,
	%since we're going to need the k_quad points anyway, so there's no
	%point calculating them twice. Do it manually
	k_quad = unwarp_data('linear', pts, ka, kb);
	Abar = interp1(ki,A,k_quad,'linear','extrap');
	
	z_quad = unwarp_data(opts.warp_method, pts, zf);
	deriv = warp_derivative(opts.warp_method, pts, zf);
	
	[zz,kk] = meshgrid(z_quad, k_quad); %z increases across row, k decreases down column
	common_term = ((alpha^2)./kk + i*(zz-z0))./kk;
	K_unweighted = 0.5*i*diag(Abar)*((1./common_term).*exp(2*i*(zz-z0).*kk-Q^2/4*common_term))*diag(deriv);
	W = diag(weights);
	Kd = K_unweighted*W;
	Kdag = ctranspose(K_unweighted)*W;
	
end

