%Iterative solver for large 2nd kind Fredholm equation in form
% $S = \eps \chi + \hat{K} \chi
% Where K is positive definite
function [chi, error, iterations] = solve_iteratively(K, S, eps, varargin)
	%Default tolerance and maximum iterations
	tol = 0.0001;
	maxIters = 1000;
	if nargin > 3
		tol = vargin{1};
	end;
	if nargin > 4
		maxIters = vargin{2};
	end;
	
	%select relaxation parameter. 0.9 is a safety factor
	%transform equation into form $(I - \lambda \hat{K})\chi = y$
	y = S/eps;
	lambda = -1/eps;
	normK = operator_norm(K);
	mu = 0.9 * -2*lambda/(normK - lambda);
	
	%error estimates: |e| <= |x[n+1] - x[n]|/(1-|(1-mu)I + mu*lambdaK|)
	[r,c] = size(K);
	new_op = (1-mu)*eye(r) + mu*lambda*K;
	error_multiplier = 1-operator_norm(new_op);
	abs_tol = tol*error_multiplier;
	
	iterations = 0;
	error = abs_tol + 1;
	chi = mu*y;
	
	while error > abs_tol && iterations < maxIters
		chi_old = chi;
		chi = mu*y + (1 − mu)*chi + mu*lambda*K*chi;
		error = norm(chi - chi_old);
		iterations = iterations + 1;
	end;
	
	error = error/error_multiplier;

end
