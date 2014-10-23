%Solve a 3D OCT problem, including the Fourier transforms
%
%Usage:
%[chi, rx, ry, rz, error] = solve_3d(f, S, S_xi, S_yi, ki, A, zf, varargin)
%
%Input Parameters
%	f: Fastcall kernel. f should be a function handle that
%	accepts arguments Qx, Qy, and returns a 1D fastcall kernel generator
%	For example, use 
%		f = @(Qx,Qy) fastcall_gauss_kernel(sqrt(Qx^2+Qy^2),alpha,z0);
%	 to set alpha and z0 on gaussian beam
%
%	S: recorded interferometric data, S(x,y,k) in rank three tensor. 
%	Major dimension is x, then y to aid fft performance. k is stored in
%	dimension with largest stride
%	
%	S_xi, S_yi: Spatial sampling points of S. S_xi and S_yi must contain
%	uniformally spaced sampling points
%
%	ki: sampling points in k domain of both S and A
%
%	A: Source spectrum envelope. Must be specified at the ki points
%
%	zf: scatterer thickness along z-axis
%
%	varargin:
%		do_fft: By default set to 1, in which case we do the 2D fourier
%		transforms. If set to 0, it will be assumed that S is already 
%		dft, and S_xi and S_yi actually contain the frequencies
%		
%		opts_1d: additional options to pass to the 1d solver. Refer to 
%		function solve_1d_opts for allowed options and parameters
%
%Output Parameters
%	chi: susceptibility, as a rank three tensor. chi(x,y,z), x is major
%	dimension, z is minor dimension
%
%	rx, ry, rz: sampling points in spatial domain for chi.
%
%	error: not currently defined or used


function [chi, rx, ry, rz, error] = solve_3d(f, S, S_xi, S_yi, ki, A, zf, varargin)
	%get additional options for the 1d solver
	if nargin > 7
		do_fft = varargin{1};
	else
		do_fft = 1;
	end;
	if nargin > 8
		opts_3d = varargin{2};
	else
		opts_3d = solve_3d_opts(); %use defaults
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
	if (opts_3d.min_feature < pi/kmax)
		wc = 2*zf*kmax; %diffraction limit for non dimensionalised
	else
		wc = 4*pi/(opts_3d.min_feature/zf);
	end;
	
	%if no discretisation level set, set based on Diffraction limit - let
	%'Nyquist' rate be 1.5 times diffraction limit. 1.5 is somewhat arbitrary magic number,
	%just chosen different from 1 to gives us some design margin
	%Note that quad points aren't equally spaced, hence 'Nyquist' rate = 1/2max(Ts)
	if opts_3d.n == 0
		Nyquist = 1.5*wc;
		max_sample_period = 0.5/Nyquist;
		opts.n = quad_points_needed(opts_3d.quad_method,max_sample_period);
	end;
	
	%Precompute the quadrature points, weights, low pass filter as a matrix.
	[pts, weights] = generate_quadrature(opts_3d.quad_method, opts_3d.n);
	%TODO: pass computed points and weights to the 1d solver so they don't
	%have to be recomputed
	if (opts_3d.do_filter)
		P = make_lpf(wc,pts,weights);
		filter = @(x) P*x;
	else
		filter = @(x) x; %identity, no filtering
	end;
	
	%Create 1d options from the 3d options
	opts_1d = solve_1d_opts(
		'solver', opts_3d.solver,...
		'n',opts_3d.n,...
		'mean_chi',opts_3d.mean_chi,...
		'quad_method', opts_3d.quad_method,...
		'max_iters',opts_3d.max_iters,...
		'tol', opts_3d.tol,...
		'basis',opts_3d.basis,...
		'min_feature',opts_3d.min_feature,...
		'filter',filter);
	
	
	%TODO: need to deal with non-uniform sampling in xy plane before fft
	%perform fft on each slice of S
	if do_fft
		Stilde = fft2(S);
		
		%assume z axis passes through middle of object
		x_width = S_xi(end) - S_xi(1);
		if (x_width > eps) %Check non-zero (2D case)
			B_x = 2*pi/x_width;
			N_x = length(S_xi);
			%spatial frequencies from 0 to B/2 are in samples 1 to N/2.
			%spatial frequencies from -B/2 to <0 are in samples N/2+1 to N
			dQx = B_x/N_x;
			Qx = [0:dQx:B_x/2, -1*B_x/2:dQx:(-1*dQx)];
		else
			Qx = [0];
		end;
		
		y_width = S_yi(end) - S_yi(1);
		if (y_width > eps) %Might be zero in 2D case
			B_y = 2*pi/y_width;
			N_y = length(S_yi);
			dQy = B_y/N_y;
			Qy = [0:dQy:B_y/2, -1*B_y/2:dQy:(-1*dQy)];
		else
			Qy = [0];
		end;
		
	else
		Stilde = S;
		Qx = S_xi;
		Qy = S_yi;
	end;
	
	%We reshape the array S now to get the z axis as the major dimension
	%This is to optimise memory access
	[ix_max, iy_max, npoints] = size(Stilde);
	%TODO: rearrange memory to get z on major dimension to improve cache access
	%Stilde = shiftdim(Stilde,1); 
	z_pts = []; %Initialization in scope
	%these loops can be parallelised
	for ix = 1:ix_max
		for iy = 1:iy_max
			f1d = f(Qx(ix), Qy(iy));
			S1d = reshape(Stilde(ix,iy,:), npoints, 1);
			[chiz, z_pts, error] = solve_1d(f1d, S1d, A,...
				ki, zf, opts_1d );
			Stilde(ix,iy,:) = chiz; %reuse memory
		end;
	end;
	
	%perform ifft on each slice. Shift back to optimise memory access pattern
	%Stilde = shiftdim(Stilde,-1);
	chi = ifft2(Stilde);
	rx = S_xi;
	ry = S_yi;
	rz = z_pts;
	
end
