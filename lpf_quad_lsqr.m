%Low pass filter for non uniformally sampled data. Use the lsqr algo
%To solve the required equations
function psi_l = lpf_quad_lsqr(psi, pts, wc, varargin)
	if nargin > 3
		weights = varargin{1};
	else
		%No weights specified, assume equally spaced points => unit weight
		weights = ones(length(pts),1);
	end;
	
	%work out nyquist rate to make sure we're ok
	Ts_ext = max(diff(pts)); %Equivalent sampling rate of input data
	if (wc >= pi/Ts_ext)
		warning('Filter cutoff frequency (%f rad/s) is greater than Nyquist frequency (%f rad/s). Possible undersampling.', wc, pi/Ts_ext);
		psi_l = psi;
	else
		sqrtW = diag(sqrt(weights));
		Ts = pi/wc;
		N = ceil(1/Ts);
		if mod(N,2) == 0
			N = N+1;
		end;
		%Calculate the A matrix
		omega = pi*((1-N):2:(N-1));
		[ww,zz] = meshgrid(omega, pts);
		A = exp(i*ww.*zz);
		psi_l = A*lsqr(sqrtW*A, sqrtW*psi);
	end;
end
