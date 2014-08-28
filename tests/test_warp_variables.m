%Boilerplate
function ret = test_warp_variables(varargin)
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

%These functions test the warp_variables function from the museum
function [status, msg] = integrate_lorentzian()
	%\int_0^\infty \frac{k}{1+z^2} dz = \frac{k\pi}{2}
	H = @(k,z) k./(1+z.^2);
	[kbar, k, zbar, z, deriv] = warp_variables(0, 1, 1, 'atan');
	k0 = 2;
	int = quad(@(zz) H(k(k0),z(zz))*deriv(zz), 0, 1);
	[status, msg] = assert_eq(int, pi);
end

function [status, msg] = integrate_lorentzian_cutoff_lin_rescale()
	H = @(k,z) k./(1+z.^2);
	[kbar, k, zbar, z, deriv] = warp_variables(0, 1, 2000, 'linear');
	k0 = 2;
	int = quad(@(zz) H(k(k0),z(zz))*deriv(zz), 0, 1);
	[status, msg] = assert_eq(int, pi, 0.01); %cutoff is not an accurate method
end

function [status, msg] = test_k_warp()
	%should get a linear transformation from [ka, kb] to [0,1]
	H = @(k,z) k./(1+z.^2);
	[kbar, k, zbar, z, deriv] = warp_variables(-5, 7, 1);
	k_i = -5:7;
	kbar_i = 0:(1/12):1;
	[status, msg] = assert_eq(kbar(k_i), kbar_i);
end

%these functions test warp_data against warp_variables from the museum
function [status, msg] = test_warp_data_linear()
	k = 10*rand(50,1);
	kmin = min(k);
	kmax = max(k);
	kbar = warp_variables(kmin,kmax,1,'linear');
	[status, msg] = assert_eq(kbar(k), warp_data('linear',k,kmin,kmax));
end

function [status, msg] = test_warp_data_atan()
	z = 17.2*rand(50,1);
	kmin = 1;
	kmax = 100;
	alpha = 4;
	[kappa_func,k_func,zeta_func] = warp_variables(kmin,kmax,alpha,'atan');
	[status, msg] = assert_eq(zeta_func(z), warp_data('atan',z,alpha));
end

%these functions test unwarp_data against warp_variables from the museum
function [status, msg] = test_unwarp_data_linear()
	kappa = rand(50,1);
	kmin = 0;
	kmax = 10;
	[kappa_func,k_func] = warp_variables(kmin,kmax,1,'linear');
	[status, msg] = assert_eq(k_func(kappa), unwarp_data('linear', kappa, kmax));
end
function [status, msg] = test_unwarp_data_atan()
	zeta = rand(50,1);
	alpha = 2;
	[kappa_func,k_func,zeta_func,z_func] = warp_variables(0,1,alpha,'atan');
	[status, msg] = assert_eq(z_func(zeta), unwarp_data('atan', zeta, alpha));
end

%these functions test warp_derivative against warp_variables from the museum
function [status, msg] = test_warp_derivative_linear()
	kmin = 1;
	kmax = 3;
	zf = 5;
	[kbar,k,zbar,z,deriv] = warp_variables(kmin, kmax, zf, 'linear');
	zeta = rand(50,1);
	[status, msg] = assert_eq(deriv(zeta), warp_derivative('linear',zeta,zf));
end
function [status, msg] = test_warp_derivative_atan()
	kmin = 0.1;
	kmax = 100;
	zf = 1.234567;
	[kbar,k,zbar,z,deriv] = warp_variables(kmin, kmax, zf, 'atan');
	zeta = rand(50,1);
	[status, msg] = assert_eq(deriv(zeta), warp_derivative('atan',zeta,zf));
end
