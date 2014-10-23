%Tests for our implementation of the LSQR algo for octave. This test
%doesn't have to be run on matlab (which has a built in lsqr implementation)
%However, the tests should pass on both, the octave implementation is 
%designed to be matlab compatible
function ret = test_lsqr(varargin)
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

function [status, msg] = test_exact()
	%A system of equations that has an exact solution
	%We will solve quite ill conditioned systems, but we'll make
	%sure the condition number is < 100
	A = rand(5,5); %quite ill conditioned
	while cond(A) > 100
		A = rand(5,5);
	end;
	x = rand(5,1);
	b = A*x;
	x_lsqr = lsqr(A,b,0.001,15);
	[status, msg] = assert_eq(x, x_lsqr,0.001);
end

function [status, msg] = test_noisy()
	%Well conditioned system with added noise
	A = diag(1:5);
	x = rand(5,1);
	b = A*x;
	b = awgn(b, 10);
	x_exp = pinv(A)*b;
	x_lsqr = lsqr(A,b,0.001);
	[status, msg] = assert_eq(x_exp,x_lsqr,0.001);
end

function [status, msg] = test_overdetermined()
	%Overdetermined linear regression
	t = (1:20)';
	y = 2*t+1;
	y = awgn(y,1);
	A = [t, ones(20,1)];
	p = pinv(A)*y;
	p_lsqr = lsqr(A, y, 0.001);
	[status, msg] = assert_eq(p,p_lsqr,0.001);
end

function [status, msg] = test_damped()
	%Solve a damped system
	A = diag(1:5);
	epsilon = 0.1;
	x = rand(5,1);
	b = inv(A)*(A'*A*x + epsilon*x);
	b_damped = [b; zeros(5,1)];
	nk = length(b);
	nz = 5;
	%This is basically a ternary statement, implemented here because
	%matlab is a terrible language
	%iif  = @(varargin) varargin{2*find([varargin{1:2:end}], 1, 'first')}();
	f = @(x,t) damp_lsqr(A, sqrt(epsilon), x, nk, nz, t);
	x_lsqr = lsqr(f,b_damped,0.001);
	[status, msg] = assert_eq(x,x_lsqr,0.001);
end
