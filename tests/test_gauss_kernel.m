%Boilerplate
function ret = test_gauss_kernel(varargin)
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

function [status, msg] = test_kernel()
	%pick some random values of alpha, Q, z0
	Q = 3*rand(1,50);
	alpha = rand(1,50);
	z0 = 10*rand(1,50);
	
	a = zeros(50,50);
	b = zeros(50,50);
	for ii = 1:50
		handle = gauss_kernel(Q(ii), alpha(ii), z0(ii));
		z = 6*(rand(1,50)-1);
		k = 6*(rand(1,50)-1);
		a(ii,:) = arrayfun(handle, k, z);
		b(ii,:) = arrayfun(@(k,z) 2*pi*i*1/4/pi/(alpha(ii)^2/k^2+i*(z-z0(ii))/k) * exp (2*i*(z-z0(ii))*k)*exp(-1*Q(ii)^2/2*(alpha(ii)^2/2/k^2+i*(z-z0(ii))/2/k)),k,z);
	end;
	%disp(a(1:3,1:3));
	%disp(b(1:3,1:3));
	[status, msg] = assert_eq(a,b);
end;
