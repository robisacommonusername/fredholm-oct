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
	
	%Correction factors for solving mean case (i.e. only solving Q=0 case)
	%Allow corrections to be specified as a row or a column vector
	[corr_r, corr_c] = size(opts.correction);
	if corr_r > corr_c
		S = S.*opts.correction;
	else
		S = S.*transpose(opts.correction); %allow complex corrections
	end;
	
	%Determine minimum feature size, or set to diffraction limit
	if (opts.min_feature < pi/kmax)
		wc = 2*zf*kmax; %diffraction limit for non dimensionalised
	else
		wc = 4*pi/(opts.min_feature/zf);
	end;
	
	%if no discretisation level set, set based on Diffraction limit - let
	%'Nyquist' rate be 1.5 times diffraction limit. 1.5 is somewhat arbitrary magic number,
	%just chosen different from 1 to gives us some design margin
	%Note that quad points aren't equally spaced, hence 'Nyquist' rate = 1/2max(Ts)
	if opts.n == 0
		Nyquist = 1.5*wc;
		max_sample_period = 0.5/Nyquist;
		opts.n = quad_points_needed(opts.quad_method,max_sample_period);
	end;
	
	%Options to pass to the fastcall function.
	fast_opts = fastcall_opts('quad_method', opts.quad_method,'n',opts.n,'low',0);
	
	%Calculate discretised operator, quad points, etc
	[Kd,Kdag,pts,weights] = f(A,ki,zf,fast_opts);
	sigma = resample_vector(S, ki, pts);
	
	%solve equations
	normK = operator_norm(Kd,Kdag,weights);
	x0 = opts.mean_chi*ones(length(sigma),1);
	switch (opts.solver)
		case 'richardson_lpf'
			%Construct low order discretisation
			[Kd_low, Kdag_low, pts_low] = f(A,ki,zf,setfield(fast_opts, 'low', 1));
			sigma_low = resample_vector(S, ki, pts_low);
			
			[chi, error, iters] = solve_iteratively_lpf(Kd, Kdag, sigma,...
				pts, weights, wc, Kd_low, Kdag_low, sigma_low,...
				solve_iteratively_opts('x0',x0,'norm_k',normK,'tol',opts.tol,...
					'max_iters',opts.max_iters));
		case 'richardson_zero'
			%Construct low order discretisation
			[Kd_low, Kdag_low, pts_low] = f(A,ki,zf,setfield(fast_opts, 'low', 1));
			sigma_low = resample_vector(S, ki, pts_low);
			
			[chi, error, iters] = solve_iteratively(Kd, Kdag, sigma, weights,...
				Kd_low,Kdag_low,sigma_low,...
				solve_iteratively_opts('x0',x0,'norm_k',normK,'tol',opts.tol,...
					'max_iters',opts.max_iters));
		case 'richardson_w2'
			%Construct low order discretisation
			[Kd_low, Kdag_low, pts_low] = f(A,ki,zf,setfield(fast_opts, 'low', 1));
			sigma_low = resample_vector(S, ki, pts_low);
			
			[chi, error, iters] = solve_iteratively_w2(Kd, Kdag, sigma,...
				pts, weights, wc, Kd_low, Kdag_low, sigma_low,...
				solve_iteratively_opts('x0',x0,'norm_k',normK,'tol',opts.tol,...
					'max_iters',opts.max_iters));
		case 'bicg_galerkin'
			[chi, error, iters] = solve_bicg_galerkin(Kd, sigma, pts, weights, wc,...
			solve_iteratively_opts('x0',x0,'norm_k',normK,'tol',opts.tol,...
				'max_iters',opts.max_iters,'basis',opts.basis));
		
		otherwise
			warning('Unrecognised solver specified. Falling back to QR');
			%Fallback/reference solver. This actually works ok on Octave,
			%because octave automatically finds least squares solutions
			%for ill-conditioned linear systems. On Matlab, it will
			%potentially fail (Matlab naively uses LU solver regardless
			%of conditioning. Refer to "Algorithms" section of
			%http://www.mathworks.com.au/help/matlab/ref/mldivide.html).
			%
			%Failure is particularly likely on Matlab with the 
			%'lcurve_lpf' regularisation, as this method will tend to 
			%select smaller regularisation parameters. Regardless of 
			%whether it works, it will be slow, because we need to 
			%compute eps*I+Kdag*K explicitly (N^3), rather than iterative
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
