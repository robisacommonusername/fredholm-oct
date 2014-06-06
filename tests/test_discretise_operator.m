%Boilerplate
function ret = test_discretise_operator(varargin)
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

function [status, msg] = test_trivial()
	H = @(k,z) exp(-(k^2+z^2));
	[pts,weights] = generate_quadrature('trivial', 2);
	[Kd, Kdag] = discretise_operator(H, pts, weights);
	Ka = 0.5*[exp(-1/8) exp(-5/8); exp(-5/8) exp(-9/8)];
	[status, msg] = assert_eq(Ka, Kd);
end;
