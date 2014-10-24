%Boilerplate
function ret = test_solve_1d(varargin)
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

function [status, msg] = small_test_q0_gaussianbeam_10db()
	lambda = 1;
	l_min = lambda-0.5;
	l_max = lambda+0.5;
	kmin = 2*pi/l_max;
	kmax = 2*pi/l_min;
	NA = 0.2;
	alpha = pi/NA;
	fwhm = 2*sqrt(3)*alpha^2/(kmin+kmax)*2;
	Q = 0;
	zf = 900*lambda;
	z0 = zf/2;
	npoints = 800;
	ki = transpose(linspace(kmin,kmax,npoints));
	A = ones(npoints,1);
	A = A/norm(A);
	noise_ratio = 0.1; %10dB SNR
	SNR = -10*log10(noise_ratio);
	solver_opts = solve_1d_opts('mean_chi',0.825,'n',npoints,'solver','bicg_galerkin','min_feature',zf/2);
	
	%thin object, one fwhm in width. refractive index varies between
	%n=1.4 => chi = 0.96, and n=1.3 => chi = 0.69
	%chi = @(z) (heaviside(z-1*lambda)-heaviside(z-4*lambda))*((0.96-0.69)*0.5*(cos(2*pi*8*z/zf)+1)+0.69);
	chi = @(z) ((0.96-0.69)*0.5*(cos(2*pi*2*z/zf)+1)+0.69);
	bg = @(z) 0.81*(heaviside(z-1.1*lambda)-heaviside(z-3.9*lambda));
	%bg = @(z) 0;
	
	f = fastcall_gauss_kernel_lsq(Q,alpha,z0);
	[chi_exp,z_exp] = generate_and_solve(chi,f,A,ki,zf,noise_ratio,solver_opts,sprintf('Thin object, NA=%1.1f, SNR=%ddB',NA, SNR),'z/lambda');
	
	fprintf('\n------- small_test_q0_gaussianbeam_10db--------\nUser Input required:\n\n');
	ans = input('Do the two graphs (more or less) coincide? Type Y for YES or N for NO ','s');
	fprintf('\n');
	if (ans(1) == 'y' || ans(1) == 'Y')
		status = 0;
		msg = 'PASS';
	else
		status = 1;
		msg = 'FAIL';
	end;
end
