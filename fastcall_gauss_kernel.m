%This function is a prototype for implementation in C

%This function supersedes
%warp_variables, gauss_kernel, discretise_operator

function f = fastcall_gauss_kernel(Q,alpha,z0,varargin)
	if nargin > 3
		opts = varargin{1};
	else
		opts = fastcall_opts();
	end;
	
	f = @(Abar,ka,kb,zf,low) fastcall_gauss_worker(Abar,ka,kb,zf,Q,alpha,z0,low,opts);
end

%This is the fastcall function, and it should be implemented in C
function [Kd,Kdag,pts,k,z,deriv] = fastcall_gauss_worker(Abar,ka,kb,zf,Q,alpha,z0,low,opts)
	if low
		quad_method = opts.quad_method_low;
		n = opts.n_low;
	else
		quad_method = opts.quad_method;
		n = opts.n;
	end;
	
	%model implementation
	%TODO: optimise
	H = gauss_kernel(Q, alpha, z0);
	
	%remap variables
	[kbar,k,zbar,z,deriv] = warp_variables(ka,kb,zf, opts.warp_method);
	
	[pts, weights] = generate_quadrature(quad_method, n);
	
	%This step is very slow (it was what we were trying to avoid by creating
	%the fastcall interface), so rewrite it.
	[Kd, Kdag] = discretise_operator(@(kk,zz) H(k(kk), z(zz))*deriv(zz), pts, weights, Abar);
end
