%Iterative solver for large 2nd kind Fredholm equation in form
% KdagS = eps*(x-x0) + Kdag*K*x
%
function [chi, error, iterations] = solve_iteratively(Kd, Kdag, S, eps, weights, varargin)
	%specify initial guess
	if nargin > 5
		x0 = varargin{1};
	else
		x0 = zeros(length(S),1);
	end;
	
	%allow user to specify |K| to prevent recalculating it
	if nargin > 6
		if varargin{2} ~= 0
			normK = varargin{2};
		else
			normK = operator_norm(Kd,Kdag,weights);
		end;
	else
		normK = operator_norm(Kd,Kdag,weights);
	end;
	
	%Default tolerance and maximum iterations
	tol = 0.0001;
	maxIters = 10000;
	if nargin > 7
		tol = varargin{3};
	end;
	if nargin > 8
		maxIters = varargin{4};
	end;
	
	%select relaxation parameter. 0.5 is a safety factor
	%transform equation into form $(lambdaI - \hat{K})deltax = y$
	%where chi = x0+deltax
	lambda = -1*eps;
	y = -1*(Kdag*S+eps*x0) - lambda*x0 + Kdag*(Kd*x0);
	mu = 0.5 * -2*lambda/(normK^2 - lambda);
	
	%error estimate
	%M = max([1-mu ; mu*normK^2/eps + mu - 1]);
	%error_multiplier = 1-M;
	error_multiplier = mu;
	abs_tol = tol*error_multiplier;
	
	iterations = 0;
	error = abs_tol + 1;
	deltax = (mu/lambda)*y;
	while error > abs_tol && iterations < maxIters
		deltax_old = deltax;
		deltax = (mu/lambda)*y + (1-mu)*deltax + Kdag*(Kd*((mu/lambda)*deltax));
		error = sqrt(weights' * (deltax - deltax_old).^2);
		iterations = iterations + 1;
	end;
	chi = x0 + deltax;
	
	if iterations == maxIters
		disp('WARNING: solve_iteratively returned after maximum number of iterations was exceeded');
	end;
	
	error = error/error_multiplier;

end
