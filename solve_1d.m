%Solve a 1d fredholm problem
% 
%Paramaters
%H: kernel function in variables k,z
%S: recorded interferometric data as a function of k (vector)
%S_ki: sampling points of S in k domain
%A: source spectrum as function of k (vector).
%A_ki: sampling points of A(k) in k domain
%pen_depth: approximate penetration depth of beam
%mean_chi: average susceptibility, required for regularisation 
%varargin
%n: discretisation level, 0 for num points in S
%tolerance, or 0 for default=0.001. Tolerance is RMS relative tolerance
%maximum iterations, or 0 for default=1000;
%quadrature method, default = gauss10
%regularisation method: Actually the method of selecting a regularisation parameter
%	either 'disc' or 'lcurve'
%snr: SNR of S, in dB. Ignored when using lcurve. Default value is 5dB

function [chi, z_pts, error] = solve_1d(H, S, S_ki, A, A_ki, pen_depth, mean_chi, varargin)

	%default parameters
	mean_chi = 0.8125;
	quad_method = 'gauss10';
	maxIters = 1000;
	tol = 0.001;
	reg_method = 'lcurve';
	snr_dB = 5;
	
	mand_params = 7;
	%parse varargin
	if nargin >= mand_params+1 && varargin{1} ~= 0
		n = varargin{1};
	else
		n = length(S);
	end;
	if nargin >= mand_params+2 && varargin{2} ~= 0
		tol = varargin{2};
	end;
	if nargin >= mand_params+3 && varargin{3} ~= 0
		maxIters = varargin{3};
	end;
	if  nargin >= mand_params+4
		quad_method = varargin{4};
	end;
	if  nargin >= mand_params+5
		reg_method = varargin{5};
	end;
	if  nargin >= mand_params+6
		snr_dB = varargin{6};
	end;	
	
	%Calculate quadrature points and weights
	[pts, weights] = generate_quadrature(quad_method, n);
	
	%remap variable into [0,1]x[0,1]
	[Hbar, kbar, z] = warp_variables(H, A_ki(1), A_ki(end), pen_depth);
	
	%resample S and A at the appropriate sampling points
	Sbar = discretise_function(S, pts, kbar(S_ki));
	Abar = discretise_function(A, pts, kbar(A_ki));
	
	%construct discretised operator and its adjoint
	[Kd, Kdag] = discretise_operator(Hbar, pts, weights, Abar);
	
	%Compute regularisation parameter
	normK = operator_norm(Kd,Kdag,weights);
	eps = regularise2(Hbar, Abar, pts, Sbar, pts, reg_method);
	
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
	
	x0 = mean_chi*ones(length(Sbar),1);
	%now use initial solution to guide iterative solver
	[chi, error, iters] = solve_iteratively(Kd, Kdag, Sbar, eps, weights, x0, normK, tol, maxIters);
	
	%unwarp z axis
	z_pts = z(pts);
	
end
