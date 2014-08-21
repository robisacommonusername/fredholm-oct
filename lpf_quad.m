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
		X(Nc:end-Nc+1) = 0; %this will give us Gibbs?
	
		%reconstruct at the quadrature points using trig interpolating
		%polynomial. 
		%Ensure to reconstruct with frequencies in [-Fs,Fs], not [0,2Fs]
		%READ:
		%https://en.wikipedia.org/wiki/Discrete_Fourier_transform#Trigonometric_interpolation_polynomial
		if (mod(N,2) == 0)
			%even length sequence - need to handle F=0.5 separately
			k = [0:(N/2-1), 0, (-1*(N/2-1)):-1];
			[kk,t] = meshgrid(k,pts);
			nyquist_terms = X(N/2+1)*cos(N*pi*pts); %F=0.5. N/2+1 due to matlab 1 indexing
			%We have to explicitly use transpose, not ', as ' will do a 
			%complex transpose by default for a complex vector - not what
			%we want. THIS was a very annoying bug to track down
			psi_l = psi_l + 1/N*(exp(2*pi*i*kk.*t)*transpose(X)) + nyquist_terms;
		else
			%odd length sequence
			k = [0:floor(N/2), (-1*floor(N/2)):-1];
			[kk,t] = meshgrid(k,pts);
			%As before, use transpose explicitly, don't do X'
			psi_l = psi_l + 1/N*(exp(2*pi*i*kk.*t)*transpose(X));
		end;
		%keyboard();
		%calculate error
		delta_psi = psi - psi_l;
		err = trapz(pts,abs(delta_psi));
		iters = iters+1;
	end;
	if (iters == max_iters)
		disp('WARNING: Voronoi method returned after exceeding maximum number of iterations');
	end;
end
