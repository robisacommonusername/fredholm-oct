%Boilerplate
function ret = test_solve_richardson_no_eps(varargin)
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

function [status, msg] = test_eps_zero()
	Kd = diag(1:15);
	x=10*rand(15,1);
	y = Kd*x;
	x_s = solve_richardson_no_eps(Kd, Kd', y, ones(15,1));
	[status, msg] = assert_eq(x,x_s,0.01);
end

function [status, msg] = test_overdetermined()
	%Solve a least squares problem
	x = (1:20)';
	y = 2*x+1;
	y = awgn(y,0.5);
	%minimise ||y-A p||, where p are the parameters m,c in y=mx+c
	A = [x,ones(20,1)];
	p_it = solve_richardson_no_eps(A,A',y,ones(2,1),solve_iteratively_opts('max_iters',10000,'tol',0.001));
	figure;
	hold on;
	scatter(x,y);
	plot(x,p_it(1)*x+p_it(2),'r');
	hold off;
	[status, msg] = assert_eq(p_it, A\y,0.1);
end
