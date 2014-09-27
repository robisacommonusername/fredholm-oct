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
	ct = bb+2*kk; %a common term
	
	%Wrong expression below: overflows
	%Kd = 0.5*i*diag(A)*kk.*exp(i*z0*bb).*exp(-(alpha^2)*ct./kk).*heaviside(ct-(Q^2)./(4*kk));
	%We have to be careful to avoid overflows here : ct./kk can be negative with large abs value
	%As a "typical" example, we might have ct as low as -22, with -alpha^2 ~ -250 for NA=0.2
	%Then exp(-alpha^2*ct./kk) = exp(5500) = overflow => will get Inf in results
	%Because of the heaviside function, we can actually ignore all the negative
	%ct values. However, ieee arith rules implies Inf*0 = Inf, so to avoid
	%Infs in our final answer, we need to not evaluate the overflowing expressions.
	%A simple solution is just to set the offending (i.e. negative) ct values
	% to zero. However, now we will have a problem when Q=0. We want
	%heaviside(ct-(Q^2)...) = 0 for all ct < 0. But when Q=0, we will get
	%heaviside(0) for all negative ct, which by default in matlab is 0.5
	%Thus we have to ensure we set the value of heaviside(0) to 0,
	%i.e heaviside(ct-..., 0) = 0 for all ct < 0 when negative ct zeroed out
	%(see help heaviside)
	ct(ct < 0) = 0;
	Kd = 0.5*i*diag(A)*kk.*exp(i*z0*bb).*exp(-(alpha^2)*ct./kk).*heaviside(ct-(Q^2)./(4*kk),0);
	Kdag = ctranspose(Kd);
	pts = beta;
	weights = ones(N,1);
	
end

