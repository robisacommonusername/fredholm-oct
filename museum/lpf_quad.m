%points must be in interval [0,1]
function psi_l = lpf_quad(psi, pts, wc)
	%work out nyquist rate and new sampling points. Sample at double Nyquist
	Ts_ext = max(diff(pts)); %Equivalent sampling rate of input data
	Ts = 0.5*Ts_ext; %Sampling rate of interpolated data
	even_pts = 0:Ts:1;
	N = length(even_pts);
	
	%Construct psi from samples, assuming that it is bandlimited in +-Wc
	%Use modified Voronoi-Allebach algorithm, with linear interpolator
	%instead of nearest neighbour
	%https://www.math.ucdavis.edu/~strohmer/research/sampling/irsampl.html
	iters = 0;
	max_iters = 5;
	psi_l = zeros(length(pts),1);
	delta_psi = psi;
	tol = 0.01*trapz(pts,abs(psi));
	err = tol+1;
	
	Nc = floor(N*wc*Ts/2/pi);
	if (wc >= pi/Ts_ext)
		warning('Filter cutoff frequency (%f rad/s) is greater than Nyquist frequency (%f rad/s). Possible undersampling.', wc, pi/Ts_ext);
		psi_l = psi;
	else
		while (err > tol && iters < max_iters)
			psi_interp = interp1(pts, delta_psi, even_pts, 'linear','extrap');
			%perform fft on interpolated data
			X = fft(psi_interp); 
			
			%reconstruct at the quadrature points using trig interpolating
			%polynomial. 
			%Ensure to reconstruct with frequencies in [-Fs,Fs], not [0,2Fs]
			%READ:
			%https://en.wikipedia.org/wiki/Discrete_Fourier_transform#Trigonometric_interpolation_polynomial
			%Note that we only take the frequencies in the band +-wc
			a = max([1, Nc]); %Make sure we always include the DC component
			Xf = [X(1:a), X(N-Nc+1:N)];
			k = [0:(a-1), (-1*Nc):-1];
			[kk,t] = meshgrid(k,pts);
			psi_l_old = psi_l;
			psi_l = psi_l + 1/N*(exp(2*pi*i*kk.*t)*transpose(Xf));
			
			%Compute correction term for next iteration
			delta_psi = psi - psi_l;
			
			%Are the iterates converging to the lowpass solution?
			err = trapz(pts,abs(psi_l - psi_l_old)); 
			iters = iters+1;
		end;
		if (iters == max_iters)
			warning('Voronoi method returned after exceeding maximum (%d) number of iterations', max_iters);
		end;
	end;
end
