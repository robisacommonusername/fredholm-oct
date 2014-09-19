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
	NA = 1.2;
	alpha = pi/NA;
	fwhm = 2*sqrt(3)*alpha^2/(kmin+kmax)*2;
	Q = 0;
	zf = 5*lambda;
	z0 = zf/2;
	npoints = 300;
	ki = linspace(kmin,kmax,npoints);
	A = ones(npoints,1);
	A = A/norm(A);
	noise_ratio = 1; %10dB SNR
	solver_opts = solve_1d_opts('mean_chi',0.8125,'n',npoints,'solver','richardson_lpf','correction',alpha);
	
	%thin object, one fwhm in width. refractive index varies between
	%n=1.4 => chi = 0.96, and n=1.3 => chi = 0.69
	chi = @(z) 0.135*(cos(2*8*pi*z/zf)+6.111).*heaviside(zf-z).*heaviside(z);
	
	
	f = fastcall_gauss_kernel(Q,alpha,z0);
	[chi_exp,z_exp] = generate_and_solve(chi,f,A,ki,zf,noise_ratio,solver_opts,'Plot Title','z/lambda');
	
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
