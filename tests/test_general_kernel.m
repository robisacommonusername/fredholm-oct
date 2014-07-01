%Boilerplate
function ret = test_general_kernel(varargin)
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

function [status,msg] = test_gaussian_kernel()
	%We should be able to recreate the gaussian kernel
	alpha = pi;
	z0 = 0;
	Q=[0,0];
	H1 = gauss_kernel(norm(Q),alpha,z0);
	g0 = @(q,k) exp(-1*(q*q')/2*alpha^2/k^2);
	H2 = general_kernel(g0,Q,z0);
	%pick some test points
	k = rand(1,50)+1;
	z = rand(1,50);
	a = arrayfun(H1,k,z);
	b = arrayfun(H2,k,z);
	[status,msg] = assert_eq(a,b,0.001*50,'abs');
end
