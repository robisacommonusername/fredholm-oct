%points must be in interval [0,1]
function psi_l = lpf_quad(psi, pts, wc)
	%work out nyquist rate and new sampling points. Sample at double Nyquist
	Ts = 0.5*max(diff(pts));
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
	while (err > tol && iters < max_iters)
		psi_interp = interp1(pts, delta_psi, even_pts, 'linear','extrap');
		%perform fft on interpolated data
		X = fft(psi_interp); 
		%filter data in freq domain - cutoff at wc
		Nc = ceil(N*wc*Ts/2/pi); %todo - windowing, etc
		
		%reconstruct at the quadrature points using trig interpolating
		%polynomial. 
		%Ensure to reconstruct with frequencies in [-Fs,Fs], not [0,2Fs]
		%READ:
		%https://en.wikipedia.org/wiki/Discrete_Fourier_transform#Trigonometric_interpolation_polynomial
		%Note that we only take the frequencies in the band +-wc
		Xf = [X(1:Nc-1), X(N-Nc+2:N)];
		k = [0:(Nc-2), (1-Nc):-1];
		[kk,t] = meshgrid(k,pts);
		psi_l_old = psi_l;
		psi_l = psi_l + 1/N*(exp(2*pi*i*kk.*t)*transpose(Xf));
		
		
		%keyboard();
		%Compute correction term for next iteration
		delta_psi = psi - psi_l;
		
		%Are the iterates converging to the lowpass solution?
		err = trapz(pts,abs(psi_l - psi_l_old)); 
		iters = iters+1;
	end;
	if (iters == max_iters)
		disp('WARNING: Voronoi method returned after exceeding maximum number of iterations');
	end;
end
