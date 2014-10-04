%Iterative solver for large 2nd kind Fredholm equation in form
% S = K*x
% Computes regularisation parameter epsilon, and then solves
% KdagS = epsilon*(x-x0) + Kdag*K*x
%
function [chi, error, iterations] = solve_cgls(Kd, Kdag, S, weights, epsilon, varargin)
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
	
	y = Kdag*S+epsilon*x0;
	[chi, flag, error, iterations] = cgs(@(x) Kdag*(Kd*x) + epsilon*x, y, opts.tol, opts.max_iters);
	
	if flag > 0
		warning('Conjugate gradient method did not converge');
	end;

end
