%Boilerplate
function ret = test_solve_3d(varargin)
	test_names = find_tests(mfilename());
	tests={};
	for ii = 1:length(test_names) 
		tests{ii} = str2func(test_names{ii}); 
	end;
	ret = run_tests(tests, test_names);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TEST FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [status, msg] = test_cylinder()
	%In this test we simulate the cylinder problem using the simplified
	%scattering model. Later in the report, we will try to solve this 
	%problem using a more rigorous FDTD simulation (and many many GBs of
	%data)
	%Simulation parmeters
	%Beam:
	%
	%   Beam type: Gaussian
	%   Beam focal plane: Focal plane somewhere inside object is ok. At the back of the object is fine.
	%   Spectral envelope: Not important. Flat is best.
	%   Number of wavelengths sampled in spectrum: 1000. Less is possible. This is limited by number of quadrature points required to perform numerical integrations accurately.
	%   Centre wavelength: 800nm
	%	Bandwidth: +- 100nm
	%   Waist Diameter: 10 wavelengths
	%   NA: 0.2
	%	Polarization: circular (left or right). Not relevant for this simulation.
	%Scatterer:
	%A dielectric cylinder, with cylinder axis parallel to the beam axis.
	%
	%  	n: 1.3
	%  	Cylinder height: 5 wavelengths
	%  	Cylinder radius: 20 wavelengths
	
	%Simulation properties:
	%   Cells/Wavelength: 20
	%   Memory consumption: ~4GiB
	%   Spectral "energy" inclusion: 95%
	%   Number of axis positions to simulate: ~20
	%   Spatial sampling period (axis position spacing): 0.0406436 *20 wavelengths
	
	lambda = 1; %non-dimensionalise everything please. Do all calcs in terms of Lambda
	lambda_min = lambda-0.125*lambda;
	lambda_max = lambda+0.125*lambda;
	k = linspace(2*pi/lambda_max,2*pi/lambda_min,1000); %wavenumbers in spectrum
	A = ones(length(k),1); %Spectral envelope
	n = 1.3; %Refractive index
	zf = 7*lambda; %thickness of simulation region. 5lambda cylinder, +1 wavelength padding on each end
	r = 20*lambda; %cylinder radius
	n_quad = 200;
	quad_method = 'gauss10';
	NA = 0.2; %numerical apperture
	alpha = pi/NA; %w0 = alpha/k
	w0_max = alpha*lambda_max/2/pi;
	z0 = 6*lambda; %focal plane at back of object

	%Simulate axis positions to edge+2 waists
	max_radius = (r+2*w0_max);
	Ts = 0.0406436*r; %magic number for 95% inclusion
	dQ = 2*pi/Ts;
	N = ceil(2*max_radius/Ts);
	%Always use odd length sequences, easier to deal with
	if (mod(N,2) == 0)
		N = N+1;
		max_radius = max_radius+Ts/2;
	end;
	%Do the transverse fourier transform analytically. Factor of N is to
	%Correct for the 1/N that will occur when we do the ifft
	chi_tilde = @(Q,z) (n^2-1)*N*besselj(1,r*Q)/(r*Q)*(heavisides(z-lambda)-heaviside(z - 6*lambda));
	
	Qs = 1/N*dQ*[0:((N-1)/2), ((1-N)/2):-1];
	S_tilde = zeros(n_quad,N,N); %z axis is major dimension will rotate dimensions later
	fast_opts = fastcall_opts('n',n_quad,'quad_method',quad_method);
	for ix = 1:N
		for iy = 1:N
			Qx = Qs(ix);
			Qy = Qs(iy);
			Q = sqrt(Qx^2+Qy^2);
			f = fastcall_gauss_kernel(Q,alpha,z0);
			[Kd,Kdag,quad_pts] = f(A,ki,zf,fast_opts);
			chi_z = arrayfun(chi_tilde,Q,quad_pts);
			S_tilde(:,ix,iy) = Kd*chi_z;
		end;
	end;
	Stilde = shiftdim(Stilde,-1);
	
	%Now attempt to solve problem. We can set do_fft to 0
	solve_opts = solve_1d_opts('n',n_quad,'quad_method',quad_method,'mean_chi',n^2-1);
	[chi_exp, rx, ry, rz] = solve_3d(f, Stilde, Qs, Qs, k, A, zf, 0, solve_opts);
	
	%Save data in case we want it for later
	d = datestr(date(),'yyyymmdd');
	save('-hdf5',sprintf('./precalc/cylinder_solved_%s.hdf5',d), 'k','A','n','zf',...
		'r','n_quad','quad_method','NA','z0','Ts',...
		'Qs','Stilde','chi_exp','rx','ry','rz');
	
	%Plot some cross sections:
	colormap('gray');
	
	%plane x=0. Should see a rectangle
	figure;
	%Find plane in data closest to x=0
	indexes_x = 1:length(rx);
	idx_x0 = interp1(rx,indexes_x,0,'nearest');
	cross_section = chi_exp(idx_x0,:,:);
	imagesc(cross_section);
	title('Cross-section, x=0');
	colorbar();
	xlabel('y');
	ylabel('z');
	print(sprintf('./img/cylinder_x0_%s.png',d));
	print(sprintf('./img/cylinder_x0_%s.eps',d));
	
	
	%plane z = 3.5*lambda (ie circular cross section through middle of cylinder)
	figure;
	indexes_z = 1:length(rz);
	idx_z35 = interp1(rz,indexes_z,3.5*lambda,'nearest');
	cross_section = chi_exp(:,:,idx_z35);
	imagesc(cross_section);
	title('Cross-section, z = 3.5L');
	colorbar();
	xlabel('x');
	ylabel('y');
	print(sprintf('./img/cylinder_z35_%s.png',d));
	print(sprintf('./img/cylinder_z35_%s.eps',d));
	
	%plane z=6*lambda (i.e. back of cylinder)
	figure;
	idx_z6 = interp1(rz,indexes_z,6*lambda,'nearest');
	cross_section = chi_exp(:,:,idx_z6);
	imagesc(cross_section);
	title('Cross-section, z = 6L');
	colorbar();
	xlabel('x');
	ylabel('y');
	print(sprintf('./img/cylinder_z6_%s.png',d));
	print(sprintf('./img/cylinder_z6_%s.eps',d));
	
	%plane z = 0.5 (should get 0)
	figure;
	idx_z05 = interp1(rz,indexes_z,0.5*lambda,'nearest');
	cross_section = chi_exp(:,:,idx_z05);
	imagesc(cross_section);
	title('Cross-section, z = 0.5L');
	colorbar();
	xlabel('x');
	ylabel('y');
	print(sprintf('./img/cylinder_z05_%s.png',d));
	print(sprintf('./img/cylinder_z05_%s.eps',d));
	
	status = 0;
	msg = 'Complete';
	
end
