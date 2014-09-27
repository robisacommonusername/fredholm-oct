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
%	be set to 0. The snr in db is -10 log10(noise_ratio)

function [Sexp, ki] = generate_test_data(chi, f, A, ki, zf, noise_ratio, opts)

	%determine max and min wavenumber - swap if necessary
	kmin = ki(1);
	kmax = ki(end);
	if kmax < kmin
		ki = flipdim(ki);
		S = flipdim(S);
		A = flipdim(A);
		kmin = ki(1);
		kmax = ki(end);
	end;
	
	%Determine minimum feature size, or set to diffraction limit
	if (opts.min_feature < pi/kmax)
		wc = 2*zf*kmax; %diffraction limit for non dimensionalised
	else
		wc = 4*pi/(opts.min_feature/zf);
	end;
	

	fast_opts = fastcall_opts('n',opts.n,'quad_method',opts.quad_method, 'wc', wc);
	[Kd,Kdag,pts,weights] = f(A,ki,zf,fast_opts);

	switch (opts.basis)
		case 'spatial'
		psi = arrayfun(chi,unwarp_data(fast_opts.warp_method, pts, zf));
		case 'frequency'
		%Need psi in the frequency domain.
		N = length(pts);
		z = linspace(0,1,N);
		psi_z = arrayfun(chi, z);
		psi = transpose(fft(psi_z));
		otherwise
		error('Unknown kernel basis %s',opts.basis);
	end;
	
	
	sigma = Kd*psi;
	%add white noise
	if noise_ratio
		snr_db = -10*log10(noise_ratio);
		if opts.downsample
			energy = 10*log10(trapz(unwarp_data('linear',pts,kmin,kmax), abs(sigma).^2));
		else
			energy = 10*log10(trapz(ki, abs(sigma).^2));
		end;
		Sexp = awgn(sigma, snr_db, energy);
	else
		Sexp = sigma;
	end;
end
