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
		opts_1d = varargin{2};
	else
		opts_1d = solve_1d_opts(); %use defaults
	end;
	
	%TODO: need to deal with non-uniform sampling in xy plane before fft
	%perform fft on each slice of S
	if do_fft
		Stilde = fft2(S);
		
		%assume z axis passes through middle of object
		x_width = S_xi(end) - S_xi(1);
		y_width = S_yi(end) - S_yi(1);
		B_x = 2*pi/x_width;
		B_y = 2*pi/y_width;
		N_x = length(S_xi);
		N_y = length(S_yi);
		%spatial frequencies from 0 to B/2 are in samples 1 to N/2.
		%spatial frequencies from -B/2 to <0 are in samples N/2+1 to N
		dQx = B_x/N_x;
		dQy = B_y/N_y;
		Qx = [0:dQx:B_x/2, -1*B_x/2:dQx:(-1*dQx)];
		Qy = [0:dQy:B_y/2, -1*B_y/2:dQy:(-1*dQy)];
	else
		Stilde = S;
		Qx = S_xi;
		Qy = S_yi;
	end;
	
	%We reshape the array S now to get the z axis as the major dimension
	%This is to optimise memory access
	Stilde = shiftdim(Stilde,1);

	[ix_max, iy_max, npoints] = size(Stilde);
	z_pts = [];
	%these loops can be parallelised
	for ix = 1:ix_max
		for iy = 1:iy_max
			f1d = f(Qx(ix), Qy(iy));
			[chiz, z_pts, error] = solve_1d(f1d, Stilde(:,ix,iy), A,...
				ki, zf, opts_1d );
			Stilde(ix,iy,:) = chiz; %reuse memory
		end;
	end;
	
	%perform ifft on each slice. Shift back to optimise memory access pattern
	Stilde = shiftdim(Stilde,-1);
	chi = ifft2(Stilde);
	rx = S_xi;
	ry = S_yi;
	rz = z_pts;
	
end
