%Essentially same as solve_iteratively, but more useful when we are
%searching for low pass solutions. Performs a low pass filtering
%on every iteration, to constrain iterates to space of low pass solutions
%
%Iterative solver for large 2nd kind Fredholm equation in form
% KdagS = epsilon*(x-x0) + Kdag*K*x
%
%A technical note: the second kind equation has a unique solution in L^2 
%(where the solution exists). i.e choosing the regularisation parameter selects
%a solution to the first kind equation. However, if epsilon is small, then the
%second kind equation will be poorly conditioned. This solver attempts to
%help with the conditioning, by restricting iterates to w^2(B). However,
%the solution may not lie within w^2(B) for a given epsilonilon. In this case
%the solver will converge to the orthogonal projection of the L^2 solution
%onto w^2(B). We hope that this solution is 'close' to the true solution,
% i.e. ||x - Px||<<||x||
function [chi, error, iterations] = solve_cgls_lpf(Kd, Kdag, S,...
	pts, weights, epsilon, wc, varargin)
	%set up options
	if nargin > 7
		opts = varargin{1};
	else
		opts = solve_iteratively_opts(); %defaults
	end;
	x0 = opts.x0;
    if x0 == 0
		[r,c] = size(Kd);
        x0 = zeros(c,1);
    end;
	
	y = lpf_quad(Kdag*S+epsilon*x0, pts, wc);
	[chi, flag, error, iterations] = cgs(@(x) lpf_quad(Kdag*(Kd*x) + epsilon*x, pts, wc), y, opts.tol, opts.max_iters);
	
	if flag > 0
		warning('Conjugate gradient method did not converge');
	end;

end
