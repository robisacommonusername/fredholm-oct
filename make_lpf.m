%Create a matrix to perform low pass filtering on a non-uniformally
%sampled signal using pseudo-inverse
function P = make_lpf(wc, pts, varargin)
	if nargin > 2
		weights = varargin{1};
	else
		%assume equally spaced points => unit weights
		weights = ones(length(pts),1);
	end;
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
	P = A*pinv(sqrtW*A)*sqrtW;
end
