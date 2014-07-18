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
