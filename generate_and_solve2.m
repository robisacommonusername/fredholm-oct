%Generate and solve a two dimensional image problem. This function can
%Be used as a fast test of solve_3d

%We generate some "point" scatterers, do forward problem, then solve
%And create an image

%f is a fastcall generator generator
function [chi_exact, chi_exp, y, z] = generate_and_solve2(f,A,ki,zf,width,n_lines,chi_fun,varargin)
	if nargin > 7
		opts = varargin{1};
	else
		opts = solve_1d_opts('n',300,'mean_chi',0);
	end;
	
	do_save = 0;
	if nargin > 8
		fn = varargin{2};
		do_save = 1;
	end;
	
	kmin = min(ki);
	kmax = max(ki);
	nk = length(ki);
	
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
	
	%Generate actual image - make some small circles
	[pts, weights] = generate_quadrature(opts.quad_method, opts.n);
	n_z = length(pts);
	z = linspace(0,zf,n_z);
	%z = unwarp_data('linear', pts, zf);
	[yy,zz] = meshgrid(y,z);
	chi_exact = chi_fun(zz,yy);
	figure;
	colormap('gray');
	imagesc(yy,zz,chi_exact);
	title('Exact Scatterer Distribution');
	xlabel('y/lambda');
	ylabel('z/lambda');
	if do_save
		print(sprintf('%s_true.png',fn));
		print(sprintf('%s_true.pdf',fn));
	end;
	keyboard();
	%Now "OCT" it
	fast_opts = fastcall_opts('n',opts.n,'quad_method',opts.quad_method);
	lines_fft = fft(chi_exact,n_lines,2);
	oct_img_fft = zeros(nk,n_lines);
	for line = 1:n_lines
		f1d = f(0,Qy(line));
		[Kd, Kdag] = f1d(A,ki,zf,fast_opts);
		oct_img_fft(:,line) = Kd*lines_fft(:,line);
	end;
	
	oct_img = ifft(oct_img_fft,n_lines,2);
	%Invert - need to have x,y on major dimensions, k on minor dimension
	S = reshape(transpose(oct_img), 1, n_lines, nk);
	[chi3d, rx, ry, rz] = solve_3d(f, S, 0, y, ki, A, zf, 1, opts);
	%Reshape chi_3d
	chi_exp = transpose(reshape(chi3d,n_lines, n_z));
	figure;
	colormap('gray');
	imagesc(yy,zz,abs(chi_exp(:,2:n_lines)));
	title('Recovered Scatterer Distribution');
	xlabel('y/lambda');
	ylabel('z/lambda');
	if do_save
		print(sprintf('%s_recovered.png',fn));
		print(sprintf('%s_recovered.pdf',fn));
	end;
	keyboard();
	
end
