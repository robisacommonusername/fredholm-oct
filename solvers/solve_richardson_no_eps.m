%Iterative solver for large equation in form
% Kdag*S = Kdag*K*x

%
function [x, error, iterations] = solve_richardson_no_eps(Kd, Kdag, S, weights, varargin)
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
	%where x = x0+deltax
	omega = 1/normK^2;
	y = Kdag*S;
	
	%error estimate
	error_multiplier = 1/(1-abs(1-omega*normK^2));
	abs_tol = opts.tol*error_multiplier;
	
	iterations = 0;
	error = abs_tol + 1;
	x = omega*y;
	while error > abs_tol && iterations < opts.max_iters
		x_old = x;
		x = omega*y + x - Kdag*(Kd*(omega*x));
		error = sqrt(weights' * (x - x_old).^2);
		iterations = iterations + 1;
	end;
	
	if iterations == opts.max_iters
		warning('solve_richardson_no_eps returned after maximum number of iterations was exceeded');
	end;
	
	error = error/error_multiplier;

end
