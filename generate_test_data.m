%Create dummy interferometric data in 1 dimension
%DEPRECTATE THIS FUNCTION??
%
%params
%Return values
%Sexp: "experimentally" generated interferometric data, as a function of k (fixed Q)
%k_i: sampling points in k domain
%
%inputs
%psi: function handle describing non-dimensionalised susceptibility on [0,1]
%f: fastcall kernel function
%A: spectrum envelope, sampled in k domain
%ki: sampling points in k domain
%zf: sample thickness
%noise_ratio: noise power/signal power. Specify this way to allow noise to
%	be set to 0. The snr in db is 10 log10(noise_ratio)
%
%use simpson quadrature, as this gives us equally spaced sampling points (more realistic)

function [Sexp, ki] = generate_test_data(chi, f, A, ki, zf, noise_ratio, opts)

	fast_opts = fastcall_opts('n',opts.n,'quad_method',opts.quad_method);
	[Kd,Kdag,pts,kfunc,zfunc,deriv] = f(A,ki,zf,fast_opts);
	ki = kfunc(pts);
	psi = arrayfun(chi,zfunc(pts));
	
	sigma = Kd*psi;
	kmin = min(ki);
	kmax = max(ki);
	
	%add white noise
	if noise_ratio
		snr_db = -10*log10(noise_ratio);
		energy = 10*log10((1/(kmax-kmin))*weights' * abs(sigma).^2);
		Sexp = awgn(sigma, snr_db, energy);
	else
		Sexp = sigma;
	end;
end
