%Essentially same as solve_iteratively, but more useful when we are
%searching for low pass solutions. Solves the w^2 regularised equation
%
%Iterative solver for large 2nd kind Fredholm equation in form
% KdagS = epsilon*x + [gamma*(I-P) + Kdag*K]x
%
%A technical note: the second kind equation has a unique solution in L^2 
%(where the solution exists). i.e choosing the regularisation parameter selects
%a solution to the first kind equation. However, if epsilon is small, then the
%second kind equation will be poorly conditioned. This solver attempts to
%help with the conditioning, by restricting iterates to w^2(B). However,
%the solution may not lie within w^2(B) for a given echilon. In this case
%the solver will converge to the orthogonal projection of the L^2 solution
%onto w^2(B). We hope that this solution is 'close' to the true solution,
% i.e. ||x - Px||<<||x||
function [chi, error, iterations] = solve_iteratively_w2(Kd, Kdag, S,...
	pts, weights, epsilon, wc, gamma,varargin)
	%set up options
	if nargin > 8
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
	
	%select relaxation parameter
	mu = epsilon/(normK^2+gamma+epsilon);
	y = (mu/epsilon) * (Kdag*S);
	
	error_multiplier = mu;
	abs_tol = opts.tol*error_multiplier;
	
	iterations = 0;
	error = abs_tol + 1;
	chi = x0;
	while error > abs_tol && iterations < opts.max_iters
		chi_old = chi;
		chi = y + (1-mu)*chi - Kdag*(Kd*((mu/epsilon)*chi)) - mu*gamma/epsilon*chi + mu*gamma/epsilon*lpf_quad_lsqr(chi, pts, wc, weights);
		error = sqrt(weights' * abs(chi - chi_old).^2);
		iterations = iterations + 1;
	end;
	
	if iterations == opts.max_iters
		warning('solve_iteratively returned after maximum number of iterations was exceeded');
	end;
	
	error = error/error_multiplier;

end
