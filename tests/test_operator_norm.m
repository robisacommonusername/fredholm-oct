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

function [status,msg] = test_real_L2_equiv()
	%If we use [1;1;...] as the quadrature weights, we whould compute
	%the same result as returned by norm(A)
	results_testfn = zeros(50,1);
	results_matlab = zeros(50,1);
	for ii = 1:50
		A = rand(5,5);
		results_testfn(ii) = operator_norm(A,A',[1;1;1;1;1]);
		results_matlab(ii) = norm(A);
	end;
	[status, msg] = assert_eq(results_testfn, results_matlab);
end;

function [status, msg] = test_complex_L2_equiv()
	%same as previous test, but use complex matrices
	results_testfn = zeros(50,1);
	results_matlab = zeros(50,1);
	for ii = 1:50
		A = rand(5,5) + i*rand(5,5);
		results_testfn(ii) = operator_norm(A,ctranspose(A),[1;1;1;1;1]);
		results_matlab(ii) = norm(A);
	end;
	[status, msg] = assert_eq(results_testfn, results_matlab);
end;

function [status, msg] = test_alternative_quadratures()
	%we should be able to use alternative quadratures, but get similar
	%results for the operator norm
	H = @(k,z) abs(k-z);
	[Kd1, Kdag1,pts1, weights1] = discretise_operator(H, 100, 'gauss10');
	[Kd2, Kdag2,pts2, weights2] = discretise_operator(H, 151, 'simpson');
	norm_gauss = operator_norm(Kd1,Kdag1,weights1);
	norm_simps = operator_norm(Kd2,Kdag2,weights2);
	[status, msg] = assert_eq(norm_gauss, norm_simps);
end;

function [status, msg] = test_convergence()
	%ensure convergence when dominant eigenvalues are very close
	tol = 0.001;
	maxIters = 1000;
	for eps = [1 0.1 0.01 0.001]
		A = diag([5+eps, 5-eps, 0.4]);
		[norm, iters] = operator_norm(A,A,[1;1;1],tol,maxIters);
		if iters >= maxIters
			status = 1;
			msg = sprintf('FAIL - did not converge for eps=%f', eps);
			return;
		end;
	end;
	status = 0;
	msg = 'PASS';
	
end;
