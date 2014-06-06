%Iterative solver for large 2nd kind Fredholm equation in form
% KdagS = eps*x + Kdag*K*x
% Where K is positive definite
function [chi, error, iterations] = solve_iteratively(K, Kdag, S, eps, weights, varargin)
	%Default tolerance and maximum iterations
	tol = 0.0001;
	maxIters = 1000;
	if nargin > 5
		tol = vargin{1};
	end;
	if nargin > 6
		maxIters = vargin{2};
	end;
	
	%select relaxation parameter. 0.9 is a safety factor
	%transform equation into form $(I - \lambda \hat{K})\chi = y$
	y = S/eps;
	lambda = -1/eps;
	normK = operator_norm(K, weights);
	mu = 0.9 * -2*lambda/(normK - lambda);
	
	%error estimates: |e| <= |x[n+1] - x[n]|/(1-|(1-mu)I + mu*lambdaK|)
	[r,c] = size(K);
	new_op = (1-mu)*eye(r) + mu*lambda*K;
	error_multiplier = 1-operator_norm(new_op, weights);
	abs_tol = tol*error_multiplier;
	
	iterations = 0;
	error = abs_tol + 1;
	chi = mu*y;
	
	while error > abs_tol && iterations < maxIters
		chi_old = chi;
		chi = mu*y + (1-mu)*chi + mu*lambda*(K*chi);
		error = sqrt(weights' * (chi - chi_old).^2);
		iterations = iterations + 1;
	end;
	
	error = error/error_multiplier;

end
