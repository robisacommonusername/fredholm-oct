%Solve a 1d fredholm problem, using fastcall object
% 
%Paramaters
%f: fastcall function
%S: recorded interferometric data as a function of k (vector)
%S_ki: sampling points of S in k domain
%A: source spectrum as function of k (vector).
%A_ki: sampling points of A(k) in k domain
%pen_depth: approximate penetration depth of beam
%mean_chi: average susceptibility, required for regularisation 
%varargin
%create a solve_1d_opts structure to pass other params to the solver.
% refer to solve_1d_opts.m for allowed fields

function [chi, z_pts, error] = solve_1d(f, S, S_ki, A, A_ki, pen_depth, varargin)
	%set up solver options
	if nargin > 6
		opts = varargin{1};
	else
		opts = solve_1d_opts();
	end;
	%if no discretisation level set, use number of points in S
	if opts.n == 0
		opts.n = length(S);
	end;
	
	%Calculate quadrature points and weights
	[pts, weights] = generate_quadrature(opts.quad_method, opts.n);
	
	%remap variable into [0,1]x[0,1]
	ka = A_ki(1);
	kb = A_ki(end);
	zf = pen_depth;
	[kbar,k, zbar, z, deriv] = warp_variables(ka,kb,zf);
	
	%resample S and A at the appropriate sampling points
	Sbar = discretise_function(S, pts, kbar(S_ki));
	Abar = discretise_function(A, pts, kbar(A_ki));
	
	%construct discretised operator and its adjoint
	fast_opts = fastcall_opts('quad_method', opts.quad_method,'n',opts.n);
	[Kd, Kdag] = f(Abar,ka,kb,zf,0,fast_opts);
	
	%Construct low order discretisation of K
	[pts_low,weights_low] = generate_quadrature(fast_opts.quad_method_low,fast_opts.n_low);
	Sbar_low = discretise_function(S, pts_low, kbar(S_ki));
	Abar_low = discretise_function(A, pts_low, kbar(A_ki));
	[Kd_low, Kdag_low] = f(Abar_low,ka,kb,zf,1);
	
	%Compute regularisation parameter
	normK = operator_norm(Kd,Kdag,weights);
	eps = regularise2(Kd_low, Kdag_low, Sbar_low, opts.reg_method);
	
	%solve equations
	
	%first find an approximate solution using low order discretisation
	%and gaussian elimination. Select order of discretisation
	%n_low = max([floor(n^(2/3)/10), 10]);
	%[pts_low, weights_low] = generate_quadrature(quad_method, n_low);
	%Abar_low = discretise_function(A, pts_low, kbar(A_ki));
	%[Kd_low, Kdag_low] = discretise_operator(Hbar, pts_low, weights_low, Abar_low);
	%LHS = eps*eye(n_low)+Kdag_low*Kd_low;
	%Sbar_low = discretise_function(S, pts_low, kbar(S_ki));
	%average_susceptibility = 0.825;
	%RHS = Kdag_low*Sbar_low + eps*mean_chi*ones(n_low,1);
	%LHS = eye(n_low) + Kdag_low*Kd_low;
	%x_low = LHS\RHS;
	%x0 = discretise_function(x_low, pts, pts_low);
	
	x0 = opts.mean_chi*ones(length(Sbar),1);
	%now use initial solution to guide iterative solver
	[chi, error, iters] = solve_iteratively(Kd, Kdag, Sbar, eps, weights,...
		solve_iteratively_opts('x0',x0,'norm_k',normK,'tol',opts.tol,...
			'max_iters',opts.max_iters));
	
	%unwarp z axis
	z_pts = z(pts);
	
end
