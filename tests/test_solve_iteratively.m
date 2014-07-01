%Boilerplate
function ret = test_solve_iteratively(varargin)
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

function [status,msg] = test_solver()
	x = rand(3,1);
	Kd = [1 3 2; 6 1 0; 1 0 9];
	eps = 0.2;
	weights = [1;1;1];
	Kdag = ctranspose(Kd);
	S = inv(Kdag)*(eps*eye(3)+Kdag*Kd)*x;
	x_s = solve_iteratively(Kd, Kdag, S, eps, weights);
	[status, msg] = assert_eq(x,x_s);
end

function [status,msg] = test_with_init()
	x = rand(5,1);
	x0 = awgn(x,20);
	Kd = rand(5,5);
	eps = 0.02;
	weights = ones(5,1);
	S = inv(Kd')*(eps*(x-x0) + Kd'*Kd*x);
	x_s = solve_iteratively(Kd, Kd', S, eps, weights,...
		solve_iteratively_opts('x0',x0,'tol',0.0001,'max_iters',10000));
	[status, msg] = assert_eq(x,x_s,0.0001);
end
