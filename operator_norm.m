%find the largest singular value of discretised Fredholm operator K
% and thus estimate the L2 norm

%x is defined at the quadrature points, so norm(x)
%may not be a good estimate of |x|. In order to compute |x|
%we need to know the quadrature weights, and |x| = sqrt(weights' * x.^2)
function [sqrt_sigma, numIters] = operator_norm(Kd,weights,varargin)
	%user can optionally specify a tolerance and max number of iterations,
	%otherwise use default values
	maxIters = 1000;
	tol = 0.0001;
	if nargin > 3
		tol = varargin{1};
	end;
	if nargin > 4
		maxIters = varargin{2};
	end;
	
	KdagK = ctranspose(Kd)*Kd;
	
	[rows,cols] = size(KdagK);
	%We expect that the singular function corresponding to the largest
	%singular value will be the least oscillatory of the singular
	%functions, hence we choose a constant x as our initial guess.
	x = ones(cols,1);
	x = x/sqrt(weights' * x.^2);
	eps = tol+1;
	numIters = 0;
	xold = x;
	
	while (eps > tol) && (numIters < maxIters)
		x = KdagK*x;
		x = x/sqrt(weights' * x.^2);
		eps = sqrt(weights' * (x-xold).^2);
		xold = x;
		numIters = numIters + 1;
	end;
	
	y = KdagK*x;
	%For complex matrices, there may be a small imaginary component
	%which shouldn't be there (ie the iterations will terminate before
	%the imaginary component can be reduced to zero). The call to abs
	%will ensure that we always return a real number
	sqrt_sigma = abs(weights' * y.^2)^0.25;
	
	
end
