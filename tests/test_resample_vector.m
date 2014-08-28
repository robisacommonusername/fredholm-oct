function ret = test_resample_vector(varargin)
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

function [status, msg] = test_simple_case()
	%A simple case that we'll work out by hand
	k = 1:11;
	S = k.^2;
	kmin = 1;
	kmax = 11;
	%S is defined at kappa = 0,0.1,0.2...1
	pts = 0.05:0.1:0.95;
	sigma_exact = ((10*(pts-0.05)+1).^2 + (10*(pts+0.05)+1).^2)./2;
	sigma_computed = resample_vector(S, k, pts);
	%keyboard();
	[status, msg] = assert_eq(sigma_exact, sigma_computed);
end
