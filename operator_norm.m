%find the largest singular value of discretised Fredholm operator K
% and thus estimate the l2 norm
function sqrt_sigma = operator_norm(Kd,varargin)
	%user can optionally specify a tolerance and max number of iterations,
	%otherwise use default values
	maxIters = 1000;
	eps = 0.0001;
	if nargin > 1
		eps = varargin{1};
	end;
	if nargin > 2
		maxIters = varargin{2};
	end;
	
	KdagK = ctranspose(Kd)*Kd;
	
	x = ones(columns(KdagK),1);
	x = x/norm(x);
	diff = eps+1;
	guard = 0;
	xold = x;
	
	while (diff > eps) && (guard < maxIters)
		x = KdagK*x;
		x = x/norm(x);
		diff = norm(x-xold);
		xold = x;
		guard = guard + 1;
	end;
	
	sqrt_sigma = sqrt(norm(KdagK*x));
	
	
end
