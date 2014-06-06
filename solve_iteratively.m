%Iterative solver for large 2nd kind Fredholm equation in form
% KdagS = eps*x + Kdag*K*x
% Where K is positive definite
function [chi, error, iterations] = solve_iteratively(Kd, Kdag, S, eps, weights, varargin)
	%allow user to specify |K| to prevent recalculating it
	if nargin > 5
		if varargin{1} ~= 0
			normK = varargin{1};
		else
			normK = operator_norm(Kd,Kdag,weights);
		end;
	else
		normK = operator_norm(Kd,Kdag,weights);
	end;
	
	%Default tolerance and maximum iterations
	tol = 0.0001;
	maxIters = 1000;
	if nargin > 6
		tol = vargin{2};
	end;
	if nargin > 7
		maxIters = vargin{3};
	end;
	
	%select relaxation parameter. 0.5 is a safety factor
	%transform equation into form $(I - \lambda \hat{K})\chi = y$
	y = -1*(Kdag*S);
	lambda = -1*eps;
	mu = 0.5 * -2*lambda/(normK^2 - lambda);
	
	%error estimate
	error_multiplier = mu*(1+normK^2/eps);
	abs_tol = tol*error_multiplier;
	
	iterations = 0;
	error = abs_tol + 1;
	%TODO: improve starting guess. Will converge correctly anyway, due
	%to contraction mapping theorem (there can only be one fixed point)
	chi = (mu/lambda)*y;

	while error > abs_tol && iterations < maxIters
		chi_old = chi;
		chi = (mu/lambda)*y + (1-mu)*chi + Kdag*(Kd*((mu/lambda)*chi));
		error = sqrt(weights' * (chi - chi_old).^2);
		iterations = iterations + 1;
	end;
	
	error = error/error_multiplier;

end
