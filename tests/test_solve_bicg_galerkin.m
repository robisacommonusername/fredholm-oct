%Boilerplate
function ret = test_solve_bicg_galerkin(varargin)
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

function [status, msg] = test_lowpass_sine_solution()
	%Generate a system of equations which has a slowly varying solution
	%In this test case, the filtered solution should be almost (to
	%within 1% or so) the same as the output of solve_iteratively
	
	t = (0:0.01:1)';
	x = sin(2*pi*t);
	N = length(t);
	weights = 1/N*ones(N,1);
	%This is just any old positive definite matrix, but we'll choose one
	%that we know is well conditioned
	K = diag(1:N);
	y = K*x;
	wc = 4*pi;
	xp = solve_bicg_galerkin(K, y, t, t, weights, wc);
	figure;
	hold on;
	plot(t,xp,'r');
	plot(t,x);
	hold off;
	fprintf('\n------- test_lowpass_sine_solution--------\nUser Input required:\n\n');
	ans = input('Do the two graphs (more or less) coincide? Type Y for YES or N for NO ','s');
	fprintf('\n');
	if (ans(1) == 'y' || ans(1) == 'Y')
		status = 0;
		msg = 'PASS';
	else
		status = 1;
		msg = 'FAIL';
	end;
	
end

function [status, msg] = test_lowpass_cosine_solution()
	%As above, but use a cosine instead of a sine to make sure the re/imag
	%parts are treated well independently
	
	t = (0:0.01:1)';
	x = cos(2*pi*t);
	N = length(t);
	weights = 1/N*ones(N,1);
	%This is just any old positive definite matrix, but we'll choose one
	%that we know is well conditioned
	K = diag(1:N);
	y = K*x;
	wc = 4*pi;
	xp = solve_bicg_galerkin(K, y, t, t, weights, wc);
	figure;
	hold on;
	plot(t,xp,'r');
	plot(t,x);
	hold off;
	fprintf('\n------- test_lowpass_cosine_solution--------\nUser Input required:\n\n');
	ans = input('Do the two graphs (more or less) coincide? Type Y for YES or N for NO ','s');
	fprintf('\n');
	if (ans(1) == 'y' || ans(1) == 'Y')
		status = 0;
		msg = 'PASS';
	else
		status = 1;
		msg = 'FAIL';
	end;
	
end

function [status, msg] = test_projected_solution()
	%In this case, we create a system of equations that has as its solution
	%a function with high frequency content, e.g. a sinc.
	%We then test that the solve returns something close to the low pass
	%filtered version of the solution (i.e. the projection into w^2(B) )
	%If the true solution is a sinc, then we should get a wider sinc.
	
	%A low pass filtered sinc should give us a sinc at a lower frequency
	fs = 3;
	Ts = 1/fs;
	t = (-10:Ts:10)';
	tbar = (t+10)/20;
	x = sinc(t);
	
	N = length(t);
	weights = 1/N*ones(N,1);
	%This is just any old positive definite matrix, but we'll choose one
	%that we know is well conditioned
	K = diag(1:N);
	eps = 100*N;
	y = K*x;
	wc = 10*pi;
	%Don't regularise the solution
	x_l = solve_bicg_galerkin(K, y, tbar, tbar, weights, wc, solve_iteratively_opts('regularise',0));

	%x_l should go approximately like 0.5*sinc(0.5t).
	%Appears to, except near the ends we have a bit too much attentuation
	%But very good agreement in the centre
	figure;
	hold on;
	plot(t,sinc(0.5*t),'r');
	plot(t,x_l);
	hold off;
	fprintf('\n------- test_projected_solution--------\nUser Input required:\n\n');
	ans = input('Do the two graphs (more or less) coincide? Type Y for YES or N for NO ','s');
	fprintf('\n');
	if (ans(1) == 'y' || ans(1) == 'Y')
		status = 0;
		msg = 'PASS';
	else
		status = 1;
		msg = 'FAIL';
	end;
	
end


