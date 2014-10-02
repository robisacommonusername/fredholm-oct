%Same as fastcall_gauss_kernel, but doesn't resample in k space.
%This kernel can be used to solve the overdetermined least squares
%problem

function f = fastcall_gauss_kernel_lsq_filon(Q,alpha,z0)
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
	
	[pts, weights] = filon_quadrature(n, ki);
	
	z_quad = unwarp_data(opts.warp_method, pts, zf);
	deriv = warp_derivative(opts.warp_method, pts, zf);
	
	[zz,kk] = meshgrid(z_quad, ki); %z increases across row, k decreases down column
	common_term = ((alpha^2)./kk + i*(zz-z0))./kk;
	K_unweighted = 0.5*i*diag(A)*((1./common_term).*exp(2*i*(zz-z0).*kk-Q^2/4*common_term))*diag(deriv);
	Kd = K_unweighted.*weights;
	Kdag = ctranspose(Kd); %This doesn't have a direct interpretation in
	%terms of the integral operator, it's the adjoint of the discretised
	%operator only.
	
end

