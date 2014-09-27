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
	
	%Set maximum iterations to size of matrix if it hasn't been set
	if opts.max_iters == 0
		opts.max_iters = opts.n;
	end;
	
	%Options to pass to the fastcall function.
	fast_opts = fastcall_opts('quad_method', opts.quad_method,'n',opts.n,'low',0, 'wc', wc);
	
	%Calculate discretised operator, quad points, etc
	[Kd,Kdag,pts,weights] = f(A,ki,zf,fast_opts);
	
	%Downsample S if required
	if opts.downsample
		sigma = resample_vector(S, ki, pts);
	else
		sigma = S;
	end;
	
	%solve equations
	normK = operator_norm(Kd,Kdag,weights);
	x0 = opts.mean_chi*ones(length(pts),1);
	%TODO:
	%ammend all solvers to handle overdetermined problems. Will need to
	%pass additional ki parameter. This way we can make downsampling optional
	switch (opts.solver)
		case 'richardson_lpf'
			%Construct low order discretisation and regularise
			[Kd_low, Kdag_low, pts_low] = f(A,ki,zf,setfield(fast_opts, 'low', 1));
			sigma_low = resample_vector(S, ki, pts_low);
			ws = 2*pi*fast_opts.n_low;
			epsilon = lcurve_calculate_eps_lpf(Kd_low, Kdag_low, sigma_low, wc/ws);
			
			%downsample
			sigma = resample_vector(S, ki, pts);
			[chi, error, iters] = solve_iteratively_lpf(Kd, Kdag, sigma,...
				pts, weights, epsilon, wc,...
				solve_iteratively_opts('x0',x0,'norm_k',normK,'tol',opts.tol,...
					'max_iters',opts.max_iters));
		case 'richardson_zero'
			%Construct low order discretisation and regularise
			[Kd_low, Kdag_low, pts_low] = f(A,ki,zf,setfield(fast_opts, 'low', 1));
			sigma_low = resample_vector(S, ki, pts_low);
			ws = 2*pi*fast_opts.n_low;
			epsilon = lcurve_calculate_eps(Kd_low, Kdag_low, sigma_low, wc/ws);

			[chi, error, iters] = solve_iteratively(Kd, Kdag, sigma, weights,epsilon,...
				solve_iteratively_opts('x0',x0,'norm_k',normK,'tol',opts.tol,...
					'max_iters',opts.max_iters));
		case 'richardson_w2'
			%Construct low order discretisation and regularise
			[Kd_low, Kdag_low, pts_low] = f(A,ki,zf,setfield(fast_opts, 'low', 1));
			sigma_low = resample_vector(S, ki, pts_low);
			ws = 2*pi*fast_opts.n_low;
			epsilon = lcurve_calculate_eps_lpf(Kd_low, Kdag_low, sigma_low, wc/ws);
			gamma = 0.5*normK^2;

			[chi, error, iters] = solve_iteratively_w2(Kd, Kdag, sigma,...
				pts, weights, epsilon, wc, gamma,
				solve_iteratively_opts('x0',x0,'norm_k',normK,'tol',opts.tol,...
					'max_iters',opts.max_iters));
		case 'bicg_galerkin'
			%This function will handle the computation of its own regularisation
			%no downsampling required, this function can handle overdetermined
			%problems. In future, change S to sigma, at the moment we have
			%downsampling on by default
			%spatial_pts = linspace(0,1,opts.n);
			%If the kernel is a frequency domain kernel, pts will not be meaningful
			if strcmp(opts.basis, 'frequency')
				pts = linspace(0,1,opts.n);
			end;
			[chi, error, iters] = solve_bicg_galerkin(Kd, sigma, ki, pts, weights, wc,...
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
			%TODO: update comments here, they are out of date
			%Actually calculate the regularisation parameter epsilon
			b = Kdag*sigma + epsilon*x0;
			A = epsilon*eye(opts.n) + Kdag*Kd;
			error = 0;
			iters = 1;
			chi = A\b;
	end;
	
	%Apply correction factors
	%Correction factors for solving mean case (i.e. only solving Q=0 case)
	%Allow corrections to be specified as a row or a column vector
	[corr_r, corr_c] = size(opts.correction);
	if corr_r > corr_c
		chi = x0 + opts.correction.*(chi-x0);
	else
		chi = x0 + transpose(opts.correction).*(chi-x0); %allow complex corrections
	end;
	
	%unwarp z axis
	z_pts = unwarp_data(fast_opts.warp_method, pts, zf);
end
