%Solve the damped least squares problem using Paige and Saunders LSQR
%algorithm

function [chi, error, iterations] = solve_lsqr(Kd, S, epsilon, varargin)
	%set up options
	if nargin > 3
		opts = varargin{1};
	else
		opts = solve_iteratively_opts(); %defaults
	end;
	x0 = opts.x0;
    if x0 == 0
		[r,c] = size(Kd);
        x0 = zeros(c,1);
    end;
	
	y = Kdag*S+epsilon*x0;
	y_damped = [y;sqrt(epsilon)*x0];
	f = @(x,t) damp_lsqr(Kd, sqrt(epsilon), x, t);
	[chi, flag, error, iterations] = lsqr(f, y_damped, opts.tol, opts.max_iters);
	
	if flag > 0
		warning('LSQR method did not converge');
	end;

end
