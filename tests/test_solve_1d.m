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
	kmin = 1;
	kmax = 3;
	alpha = 0.2;
	fwhm = 2*sqrt(3)*alpha^2/(kmin+kmax)*2;
	Q = 0;
	z0 = fwhm/2;
	zf = fwhm;
	npoints = 300;
	ki = kmin:((kmax-kmin)/(npoints-1)):kmax;
	A = ones(npoints,1);
	A = A/norm(A);
	noise_ratio = 0.1; %10dB SNR
	solver_opts = solve_1d_opts('mean_chi',0.8125,'n',npoints);
	
	%thin object, one fwhm in width. refractive index varies between
	%n=1.4 => chi = 0.96, and n=1.3 => chi = 0.69
	chi = @(z) 0.135*(cos(2*pi*z/fwhm)+6.111).*heaviside(fwhm-z).*heaviside(z);
	
	
	f = fastcall_gauss_kernel(Q,alpha,z0);
	generate_and_solve(chi,f,A,ki,zf,noise_ratio,solver_opts);
	
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
