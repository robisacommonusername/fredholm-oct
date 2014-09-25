%Iterative solver for large 2nd kind Fredholm equation in form
% S = K*x
% Computes regularisation parameter epsilon, and then solves
% KdagS = epsilon*(x-x0) + Kdag*K*x
%
function [chi, error, iterations] = solve_iteratively(Kd, Kdag, S, weights, epsilon, varargin)
	%set up options
	if nargin > 5
		opts = varargin{1};
	else
		opts = solve_iteratively_opts(); %defaults
	end;
	x0 = opts.x0;
    if x0 == 0
		[r,c] = size(Kd);
        x0 = zeros(c,1);
    end;

	if opts.norm_k == 0
		normK = operator_norm(Kd,Kdag,weights);
	else
		normK = opts.norm_k;
	end;
	
	%select relaxation parameter. 0.5 is a safety factor
	%transform equation into form $(lambdaI - \hat{K})deltax = y$
	%where chi = x0+deltax
	lambda = -1*epsilon;
	y = -1*(Kdag*S+epsilon*x0) - lambda*x0 + Kdag*(Kd*x0);
	mu = 0.5 * -2*lambda/(normK^2 - lambda);
	
	%error estimate
	%M = max([1-mu ; mu*normK^2/epsilon + mu - 1]);
	%error_multiplier = 1-M;
	error_multiplier = mu;
	abs_tol = opts.tol*error_multiplier;
	
	iterations = 0;
	error = abs_tol + 1;
	deltax = (mu/lambda)*y;
	while error > abs_tol && iterations < opts.max_iters
		deltax_old = deltax;
		deltax = (mu/lambda)*y + (1-mu)*deltax + Kdag*(Kd*((mu/lambda)*deltax));
		error = sqrt(weights' * (deltax - deltax_old).^2);
		iterations = iterations + 1;
	end;
	chi = x0 + deltax;
	
	if iterations == opts.max_iters
		warning('solve_iteratively returned after maximum number of iterations was exceeded');
	end;
	
	error = error/error_multiplier;

end
