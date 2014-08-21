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
	
	%thin object, one fwhm in width. refractive index varies between
	%n=1.4 => chi = 0.96, and n=1.3 => chi = 0.69
	chi = @(z) 0.135*(cos(2*pi*z/fwhm)+6.111).*heaviside(fwhm-z).*heaviside(z);
	z = -1*fwhm:(fwhm/300):2*fwhm;
	
	%GENERATE DATA
	%gauss kernel Q=0, with thin object. Focal plane at z=half power point
	f = fastcall_gauss_kernel(0,alpha,fwhm/2,fastcall_opts('n',301,'quad_method','simpson'));
	
	%flat spectrum, unit power
	A = ones(1,301);
	A = A/norm(A); %normalise
	
	zf = fwhm;
	[Kd,Kdag,pts,kfunc,zfunc,der] = f(A,kmin,kmax,zf,0);
	Sexp = Kd*arrayfun(chi,pts);
	ki = kfunc(pts);
	
	
	%ATTEMPT RECOVERY
	f2 = fastcall_gauss_kernel(0,alpha,fwhm/2,fastcall_opts('n',300,'quad_method','gauss10'));
	[chi_exp,z_exp] = solve_1d(f2, Sexp, A, ki, fwhm,...
		solve_1d_opts('mean_chi',0.8125,'n',300));
		
	%Plot true solution and recovered solution
	hold on;
	plot(z_exp,chi_exp);
	plot(z,arrayfun(chi,z),'r');
	hold off;
	
	status=0;
	msg = 'pass';
end
