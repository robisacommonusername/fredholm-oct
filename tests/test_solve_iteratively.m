%Boilerplate
function ret = test_operator_norm(varargin)
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

function [status,msg] = test_solver()
	x = rand(3,1);
	Kd = [1 3 2; 6 1 0; 1 0 9];
	eps = 0.2;
	weights = [1;1;1];
	Kdag = ctranspose(Kd);
	S = inv(Kdag)*(eps*eye(3)+Kdag*Kd)*x;
	x_s = solve_iteratively(Kd, Kdag, S, eps, weights);
	[status, msg] = assert_eq(x,x_s);
end;
