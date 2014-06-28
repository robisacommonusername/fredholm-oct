%Solve a 3D OCT problem, including the Fourier transforms
%
%Usage:
%[chi, rx, ry, rz, error] = solve_3d(S, S_xi, S_yi, S_ki, H, A, A_ki,...
%	pen_depth, varargin)
%
%Input Parameters
%	S: recorded interferometric data, S(x,y,k) in rank three tensor. 
%	Major dimension is x, then y to aid fft performance. k is stored in
%	dimension with largest stride
%	
%	S_xi, S_yi, S_ki: Sampling points of S. S_xi and S_yi must contain
%	uniformally spaced sampling points, S_ki need not
%
%	H: Fredholm kernel. H should be a function handle that
%	accepts arguments Qx, Qy, and returns a 1D kernel function of (k,z).
%	i.e if R represents reals, and C the complex numbers, then
%	H: (RxR) -> ((RxR) -> C)
%	For example, use 
%		H = @(Qx,Qy) gaussian_beam(sqrt(Qx^2+Qy^2),0.1,1);
%	 to set alpha and z0 on gaussian beam
%
%	A: Source spectrum envelope
%	A_ki: Sampling points in k domain of spectrum
%
%	pen_depth: expected beam penetration depth
%
%	varargin:
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


function [chi, rx, ry, rz, error] = solve_3d(S, S_xi, S_yi, S_ki, H, A, A_ki, pen_depth, varargin)
	%get additional options for the 1d solver
	if nargin > 8
		opts_1d = varargin{1};
	else
		opts_1d = solve_1d_opts(); %use defaults
	end;
	
	%TODO: need to deal with non-uniform sampling in xy plane before fft
	%perform fft on each slice of S
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

	[ix_max, iy_max, npoints] = size(Stilde);
	z_pts = [];
	%these loops can be parallelised
	for ix = 1:ix_max
		for iy = 1:iy_max
			H1d = H(Qx(ix), Qy(iy));
			[chiz, z_pts, error] = solve_1d(H1d, Stilde(ix,iy,:), S_ki,...
				A, A_ki, pen_depth, opts_1d );
			Stilde(ix,iy,:) = chiz; %reuse memory
		end;
	end;
	
	%perform ifft on each slice
	chi = ifft2(Stilde);
	rx = S_xi;
	ry = S_yi;
	rz = z_pts;
	
end
