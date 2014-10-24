%Boilerplate
function ret = test_general_kernel2(varargin)
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

function [status,msg] = flat_spec()
	%We can solve the case where g0 = 1, z0=0, Q = [0,0] and z0=0
	%analytically: should get H(k,0) = 2\pi
	g0 = @(q,k) 1;
	H = general_kernel2(g0,[0,0],0);
	a = H(1,0);
	[status,msg] = assert_eq(a,2*pi);
end

function [status,msg] = pi_square()
	%another analytically solvable case: g0 = 1/sqrt(|q|), Q=[0,0], z0=0
	%then H(k,0) = \pi^2
	g0 = @(q,k) 1/sqrt(norm(q));
	H = general_kernel2(g0,[0,0],0);
	a = H(1,0);
	[status,msg] = assert_eq(a,pi^2);
end
