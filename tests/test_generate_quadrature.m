%Boilerplate
function ret = test_generate_quadrature(varargin)
	test_names = find_tests(mfilename());
	tests={};
	for ii = 1:length(test_names) 
		tests(ii) = str2func(test_names(ii)); 
	end;
	ret = run_tests(tests, test_names);
end;

function [status, msg] = test_gauss10()
	[pts, weights] = generate_quadrature('gauss10',30);
	[status, msg] = assert_eq(sum(weights.*exp(pts)), 1.71828);
end;


function [status, msg] = test_gaussn()
	[pts, weights] = generate_quadrature('gaussn',30,0.00005);
	[status, msg] = assert_eq(sum(weights.*exp(pts)), 1.71828);
end;

function [status, msg] = test_simpson()
	[pts, weights] = generate_quadrature('simpson',31,0.00005);
	[status, msg] = assert_eq(sum(weights.*exp(pts)), 1.71828);
end;

function [status, msg] = test_trivial()
	[pts, weights] = generate_quadrature('trivial',51,0.00005);
	[status, msg] = assert_eq(sum(weights.*exp(pts)), 1.71828);
end;

function [status, msg] = test_gauss10_monotonic()
	[pts,weights] = generate_quadrature('gauss10',30,0.00005);
	[status, msg] = assert_eq(pts, sort(pts));
end;

function [status, msg] = test_gaussn_monotonic()
	[pts,weights] = generate_quadrature('gaussn',23);
	[status, msg] = assert_eq(pts, sort(pts));
end;

function [status, msg] = test_simpson_monotonic()
	[pts,weights] = generate_quadrature('simpson',17);
	[status, msg] = assert_eq(pts, sort(pts));
end;

function [status, msg] = test_trivial_monotonic()
	[pts,weights] = generate_quadrature('trivial',30);
	[status, msg] = assert_eq(pts, sort(pts));
end;
