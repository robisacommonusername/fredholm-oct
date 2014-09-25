%Generate and solve a two dimensional image problem. This function can
%Be used as a fast test of solve_3d

%We generate some "point" scatterers, do forward problem, then solve
%And create an image

%f is a fastcall generator generator
function [chi_exact, chi_exp, y, z] = generate_and_solve2(f,A,ki,zf,width,n_lines,n_points,varargin)
	if nargin > 7
		opts = varargin{1};
	else
		opts = solve_1d_opts('n',300);
	end;
	
	kmin = min(ki);
	kmax = max(ki);
	
	%Always use odd length sequence in y
	if (mod(n_lines, 2) == 0)
		n_lines = n_lines+1;
		warning('number of lines should be odd. n_lines increased to %d',n_lines);
	end;
	
	y = linspace(-1*width/2, width/2, n_lines);
	Tsy = width/n_lines;
	dQy = 2*pi/n_lines;
	omegaN = pi/Tsy;
	%spatial frequencies in order returned by matlab/octave fft
	Qy = [0:dQy:omegaN, -1*omegaN:dQy:(-1*dQy)];
	
	%Generate actual image
	chi_exact = zeros(opts.n,n_lines);
	rand_indices = sort(ceil(opts.n*n_lines*rand(1,n_points)));
	rand_susceptibilities = (0.96-0.69)*rand(1,n_points)+0.69;
	chi_exact(rand_indices) = rand_susceptibilities;
	figure;
	colormap('gray');
	imagesc(chi_exact);
	title('Exact scatterer distribution');
	xlabel('y');
	ylabel('z');
	
	%Now "OCT" it
	fast_opts = fastcall_opts('n',opts.n,'quad_method',opts.quad_method);
	lines_fft = fft(chi_exact,n_lines,2);
	for line = 1:n_lines
		f1d = f(0,Qy(line));
		[Kd, Kdag, z, weights] = f1d(A,ki,zf,fast_opts);
		lines_fft(:,line) = Kd*lines_fft(:,line);
	end;
	k_quad = unwarp_data('linear',z,kmin,kmax);
	A_quad = interp1(ki,A,k_quad,'linear','extrap');
	
	oct_img = ifft(lines_fft,n_lines,2);
	%Invert - need to have x,y on major dimensions, k on minor dimension
	S = reshape(transpose(oct_img), 1, n_lines, opts.n);
	[chi3d, rx, ry, rz] = solve_3d(f, S, 0, y, k_quad, A_quad, zf, 1, opts);
	%Reshape chi_3d
	chi_exp = transpose(reshape(chi3d,n_lines, opts.n));
	figure;
	colormap('gray');
	imagesc(abs(chi_exp(:,2:n_lines)));
	title('Recovered scatterer distribution');
	xlabel('y');
	ylabel('z');
	
end
