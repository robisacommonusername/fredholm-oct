%find the largest singular value of discretised Fredholm operator K
% and thus estimate the L2 norm
% 
%x is defined at the quadrature points, so norm(x)
%may not be a good estimate of |x|. In order to compute |x|
%we need to know the quadrature weights, and |x| = sqrt(weights' * abs(x).^2)
function [sqrt_sigma, x,numIters] = operator_norm(Kd,Kdag,weights,varargin)
	%user can optionally specify a tolerance and max number of iterations,
	%otherwise use default values
	maxIters = 1000;
	tol = 0.0001;
	if nargin > 4
		tol = varargin{1};
	end;
	if nargin > 5
		maxIters = varargin{2};
	end;
	
	[rows,cols] = size(Kd);
	%We expect that the singular function corresponding to the largest
	%singular value will be the least oscillatory of the singular
	%functions, hence we choose a constant x as our initial guess.
	x = ones(cols,1);
	eps = tol+1;
	numIters = 0;
	sigma = sqrt(weights' * abs(x).^2);
	x = x/sigma;
	
	while (eps > tol) && (numIters < maxIters)
		old_sigma = sigma;
		x = Kdag*(Kd*x);
		sigma = sqrt(weights' * abs(x).^2);
		x = x/sigma;
		eps = abs(sigma-old_sigma);
		numIters = numIters + 1;
	end;

	sqrt_sigma = sqrt(sigma);
	
	
end
