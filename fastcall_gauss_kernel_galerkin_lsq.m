%Overdetermined gaussian kernel in the frequency domain, as in
%Ralston's paper.

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
	Ts = unwarp_data(opts.warp_method, max(diff(pts)), zf);
	beta_max = pi/Ts;
	beta_min = -1*beta_max;
	
	beta_quad = unwarp_data(opts.warp_method, pts, 2*beta_max);
	deriv = warp_derivative(opts.warp_method, pts, 2*beta_max);
	
	[bb,kk] = meshgrid(beta_quad, ki); %z increases across row, k decreases down column
	ct = bb+2*kk;
	K_unweighted = 0.5*i*diag(A)*kk.*exp(i*z0*bb-(alpha^2)./(kk.*ct)).*heaviside(ct-(Q^2)./(4*kk))*diag(deriv);
	W = diag(weights);
	Kd = K_unweighted*W;
	Kdag = ctranspose(Kd);
	
end

