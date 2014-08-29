%Solve a 1d fredholm problem, using fastcall object
% 
%Paramaters
%f: fastcall function
%S: recorded interferometric data as a function of k (vector)
%A: source spectrum as function of k (vector).
%ki: sampled wavenumbers for S and A
%zf: sample thickness
%mean_chi: average susceptibility, required for regularisation 
%varargin
%create a solve_1d_opts structure to pass other params to the solver.
% refer to solve_1d_opts.m for allowed fields

function [chi, z_pts, error] = solve_1d(f, S, A, ki, zf, varargin)
	%set up solver options
	if nargin > 5
		opts = varargin{1};
	else
		opts = solve_1d_opts();
	end;
	
	%determine max and min wavenumber - swap if necessary
	kmin = ki(1);
	kmax = ki(end);
	if kmax < kmin
		ki = flipdim(ki);
		S = flipdim(S);
		A = flipdim(A);
		kmin = ki(1);
		kmax = ki(end);
	end;
	
	%if no discretisation level set, set based on Diffraction limit - let
	%'Nyquist' rate be 1.5 times diffraction limit. 1.5 is somewhat arbitrary magic number,
	%just chosen different from 1 to gives us some design margin
	%Note that quad points aren't equally spaced, hence 'Nyquist' rate = 1/2max(Ts)
	if opts.n == 0
		Nyquist = 1.5*2*kmax*zf;
		max_sample_period = 0.5/Nyquist;
		opts.n = quad_points_needed(opts.quad_method,max_sample_period);
	end;
	
	%Options to pass to the fastcall function.
	fast_opts = fastcall_opts('quad_method', opts.quad_method,'n',opts.n,'low',0);
	
	%Calculate discretised operator, quad points, etc
	[Kd,Kdag,pts,weights] = f(A,ki,zf,fast_opts);
	sigma = resample_vector(S, ki, pts);
	
	%Construct low order discretisation
	[Kd_low, Kdag_low, pts_low] = f(A,ki,zf,setfield(fast_opts, 'low', 1));
	sigma_low = resample_vector(S, ki, pts_low);
	
	%Compute regularisation parameter using low order discretisation
	%if we're using lcurve_lpf, compute diffraction limit and cutoff
	reg_opts = opts.reg_opts;
	if strcmp(opts.reg_method,'lcurve_lpf')
		wc = 2*kmax*zf; %diffraction limit for non-dimensionalised form
		%Work out spatial sampling for low order approximation
		fs = fast_opts.n_low; %Non dimensionalised version samples [0,1] => sampling freq is just number of samples (assume equally spaced)
		Wc = wc/fs;
		reg_opts = setfield(reg_opts,'Wc',Wc);
	end;
	eps = regularise2(Kd_low, Kdag_low, sigma_low, opts.reg_method, reg_opts);
	
	%solve equations
	normK = operator_norm(Kd,Kdag,weights);
	x0 = opts.mean_chi*ones(length(sigma),1);
	switch (opts.solver)
		case 'richardson_lpf'
			[chi, error, iters] = solve_iteratively_lpf(Kd, Kdag, sigma, eps, pts, weights, 2*kmax*zf,...
			solve_iteratively_opts('x0',x0,'norm_k',normK,'tol',opts.tol,...
				'max_iters',opts.max_iters));
		case 'richardson'
			[chi, error, iters] = solve_iteratively(Kd, Kdag, sigma, eps, weights,...
			solve_iteratively_opts('x0',x0,'norm_k',normK,'tol',opts.tol,...
				'max_iters',opts.max_iters));
		otherwise
			warning('Unrecognised solver specified. Falling back to QR');
			%Fallback/reference solver. This actually works ok on Octave,
			%because octave automatically finds least squares solutions
			%for ill-conditioned linear systems. On Matlab, it will
			%potentially fail (Matlab naively uses QR solver regardless
			%of conditioning). Failure is particularly likely on Matlab
			%with the 'lcurve_lpf' regularisation, as this method will
			%tend to select smaller regularisation parameters.
			%Regardless of whether it works, it will be slow, because
			%We need to compute eps*I+Kdag*K explicitly (N^3), rather than iterative
			%solvers which are N^2. There is also a loss of precision from doing
			%this explicitly
			
			%solve Kdag*sigma = eps(x-x0) + Kdag*K*x
			b = Kdag*sigma + eps*x0;
			A = eps*eye(opts.n) + Kdag*Kd;
			error = 0;
			iters = 1;
			chi = A\b;
	end
	
	
	%unwarp z axis
	z_pts = unwarp_data(fast_opts.warp_method, pts, zf);
	
end
