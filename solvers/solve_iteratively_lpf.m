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
function [chi, error, iterations] = solve_iteratively_lpf(Kd, Kdag, S,...
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

	if opts.norm_k == 0
		normK = operator_norm(Kd,Kdag,weights);
	else
		normK = opts.norm_k;
	end;
	
	%Check for case where regularisation constant is 0
	if abs(epsilon-eps) < 2*eps
		warning('No regularisation parameter. You should use solve_richardson_no_eps instead');
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
	deltax = (mu/lambda)*lpf_quad_lsqr(y,pts,wc,weights);
	while error > abs_tol && iterations < opts.max_iters
		deltax_old = deltax;
		deltax = (mu/lambda)*y + (1-mu)*deltax + Kdag*(Kd*((mu/lambda)*deltax));
		%Apply lpf
		deltax = lpf_quad_lsqr(deltax,pts,wc,weights);
		error = sqrt(weights' * abs(deltax - deltax_old).^2);
		iterations = iterations + 1;
	end;
	chi = x0 + deltax;
	
	if iterations == opts.max_iters
		warning('solve_iteratively returned after maximum number of iterations was exceeded');
	end;
	
	error = error/error_multiplier;

end
