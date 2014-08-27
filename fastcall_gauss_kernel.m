%This function is a prototype for implementation in C

%This function supersedes
%warp_variables, gauss_kernel, discretise_operator

function f = fastcall_gauss_kernel(Q,alpha,z0)
	f = @(A,ki,zf,opts) fastcall_gauss_worker(A,ki,zf,Q,alpha,z0,opts);
end

%This is the fastcall function, and it should be implemented in C
function [Kd,Kdag,pts,k,z,deriv] = fastcall_gauss_worker(A,ki,zf,Q,alpha,z0,opts)
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
	
	%model implementation
	%TODO: optimise
	H = gauss_kernel(Q, alpha, z0);
	
	%remap variables
	[kbar,k,zbar,z,deriv] = warp_variables(ka,kb,zf, opts.warp_method);
	
	[pts, weights] = generate_quadrature(quad_method, n);
	
	%resample A at the quadrature points.
	Abar = discretise_function(A, pts, ki);
	
	%This step is very slow (it was what we were trying to avoid by creating
	%the fastcall interface), so rewrite it.
	[Kd, Kdag] = discretise_operator(@(kk,zz) H(k(kk), z(zz))*deriv(zz), pts, weights, Abar);
end
