%Boilerplate
function ret = test_warp_variables(varargin)
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

function [status, msg] = integrate_lorentzian()
	%\int_0^\infty \frac{k}{1+z^2} dz = \frac{k\pi}{2}
	H = @(k,z) k/(1+z^2);
	[Hbar, kbar, z] = warp_variables(H, 0, 1, 1);
	k = 2;
	int = quad(@(zz) Hbar(k,zz), 0, 1);
	[status, msg] = assert_eq(int, pi);
end;

function [status, msg] = test_k_warp()
	%should get a linear transformation from [ka, kb] to [0,1]
	H = @(k,z) k/(1+z^2);
	[Hbar, kbar, z] = warp_variables(H, -5, 7, 1);
	k_i = -5:7;
	kbar_i = 0:(1/12):1;
	[status, msg] = assert_eq(kbar(k_i), kbar_i);
end;
