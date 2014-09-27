%Overdetermined gaussian kernel in the frequency domain, as in
%Ralston's paper.

function f = fastcall_gauss_kernel_galerkin_lsq(Q,alpha,z0)
	f = @(A,ki,zf,opts) fastcall_gauss_worker(A,ki,zf,Q,alpha,z0,opts);
end

%This is the fastcall function, and it should be implemented in C
function [Kd,Kdag,pts,weights] = fastcall_gauss_worker(A,ki,zf,Q,alpha,z0,opts)
	kmax = ki(end);
	kmin = ki(1);
	if kmin > kmax
		kmin = ki(end);
		kmax = ki(1);
	end;
	
	if abs(opts.wc) < eps
		%Set to diffraction limit
		opts.wc = 2*kmax*zf;
	end;
	
	%Calculate required N for given wc. Always use odd sequence
	Ts = pi/opts.wc;
	N = ceil(1/Ts);
	if (mod(N,2) == 0)
		N = N+1;
	end;
	Ts = 1/N;
	
	%The required beta values are (in order of octaves fft)
	beta = 2*pi*[0:((N-1)/2), ((1-N)/2):-1]';
	
	[bb,kk] = meshgrid(beta, ki); %beta increases across row, k decreases down column
	ct = bb+2*kk;
	K_unweighted = 0.5*i*diag(A)*kk.*exp(i*z0*bb-(alpha^2)./(kk.*ct)).*heaviside(ct-(Q^2)./(4*kk));
	W = diag(weights);
	Kd = K_unweighted*W;
	Kdag = ctranspose(Kd);
	
end

