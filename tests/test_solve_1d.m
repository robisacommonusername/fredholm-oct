%Boilerplate
function ret = test_solve_1d(varargin)
	test_names = find_tests(mfilename());
	tests={};
	for ii = 1:length(test_names) 
		tests(ii) = str2func(test_names(ii)); 
	end;
	ret = run_tests(tests, test_names);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TEST FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [status, msg] = small_test_q0_gaussianbeam_10db()
	kmin = 1;
	kmax = 3;
	alpha = 0.2;
	fwhm = 2*sqrt(3)*alpha^2/(kmin+kmax)*2;
	
	%thin object, one fwhm in width. refractive index varies between
	%n=1.4 => chi = 0.96, and n=1.3 => chi = 0.69
	chi = @(z) 0.135*(cos(2*pi*z/fwhm)+6.111).*heaviside(fwhm-z).*heaviside(z);
	z = -1*fwhm:(fwhm/300):2*fwhm;
	
	%gauss kernel Q=0, with thin object. Focal plane at z=half power point
	H = gauss_kernel(0,alpha,fwhm/2);
	
	%flat spectrum, unit power
	step = (kmax-kmin)/300;
	A_ki = kmin:step:kmax;
	A = ones(1,length(A_ki));
	A = A/trapz(A_ki, A); %normailse power
	
	%generate experimental data, 300 points, 
	[Sexp, k_i] = generate_test_data(301, chi, H, A, A_ki, 2*fwhm, 10);
	
	%attempt to recover solution
	[chi_exp,z_exp] = solve_1d(H, Sexp, k_i, A, A_ki, fwhm,...
		solve_1d_opts('mean_chi',0.8125,'n',300));
	plot(z_exp,chi_exp);
	
	status=0;
	msg = 'pass';
end;
