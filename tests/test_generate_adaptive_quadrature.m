%Boilerplate
function ret = test_generate_adaptive_quadrature(varargin)
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

function [status, msg] = test_gauss10()
	f = @(x) exp(-2*x)*cos(100*x);
	[pts,weights] = generate_adaptive_quadrature(f, 'gauss10', 0.00001, 10, 100);
	test = weights' * arrayfun(f, pts);
	actual = quad(f,0,1,0.00001);
	[status, msg] = assert_eq(test, actual,0.00001/abs(actual));
end;

function [status, msg] = test_gauss3()
	f = @(x) exp(-2*x)*cos(100*x);
	[pts,weights] = generate_adaptive_quadrature(f, 'gaussn', 0.00001, 3, 150);
	test = weights' * arrayfun(f, pts);
	actual = quad(f,0,1,0.00001);
	[status, msg] = assert_eq(test, actual, 0.00001/abs(actual));
end;

function [status, msg] = test_gauss3_monotonic()
	f = @(x) exp(-2*x)*cos(100*x);
	[pts,weights] = generate_adaptive_quadrature(f, 'gaussn', 0.00001, 3, 150);
	[status,msg] = assert_eq(pts, sort(pts));
end;
