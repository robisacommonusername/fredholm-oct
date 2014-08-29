function ret = test_quad_points_needed(varargin)
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

function [status, msg] = test_gauss10()
	[pts1, weights] = generate_quadrature('gauss10', 120);
	dm1 = max(diff([0;pts1;1]));
	[pts2, weights] = generate_quadrature('gauss10', 110);
	dm2 = max(diff([0;pts2;1]));
	dmax = (dm1+dm2)/2;
	n = quad_points_needed('gauss10',dmax);
	[status, msg] = assert_eq(n,120);
	
end

function [status, msg] = test_gaussn_odd()
	[pts1, weights] = generate_quadrature('gaussn', 17);
	dm1 = max(diff([0;pts1;1]));
	[pts2, weights] = generate_quadrature('gaussn', 16);
	dm2 = max(diff([0;pts2;1]));
	dmax = (dm1+dm2)/2;
	n = quad_points_needed('gaussn',dmax);
	[status, msg] = assert_eq(n,17);
end

function [status, msg] = test_gaussn_even()
	[pts1, weights] = generate_quadrature('gaussn', 18);
	dm1 = max(diff([0;pts1;1]));
	[pts2, weights] = generate_quadrature('gaussn', 17);
	dm2 = max(diff([0;pts2;1]));
	dmax = (dm1+dm2)/2;
	n = quad_points_needed('gaussn',dmax);
	[status, msg] = assert_eq(n,18);
end

function [status, msg] = test_gaussn_one()
	n = quad_points_needed('gaussn',0.75);
	[status, msg] = assert_eq(n,1);
end

function [status, msg] = test_gaussn_power2()
	[pts1, weights] = generate_quadrature('gaussn', 32);
	dm1 = max(diff([0;pts1;1]));
	[pts2, weights] = generate_quadrature('gaussn', 31);
	dm2 = max(diff([0;pts2;1]));
	dmax = (dm1+dm2)/2;
	n = quad_points_needed('gaussn',dmax);
	[status, msg] = assert_eq(n,32);
end

function [status, msg] = test_simpson_mod1()
	[pts1, weights] = generate_quadrature('simpson', 121);
	dm1 = max(diff([0;pts1;1]));
	[pts2, weights] = generate_quadrature('simpson', 119);
	dm2 = max(diff([0;pts2;1]));
	dmax = (dm1+dm2)/2;
	n = quad_points_needed('simpson',dmax);
	%keyboard();
	[status, msg] = assert_eq(n,121);
end

function [status, msg] = test_simpson_mod3()
	[pts1, weights] = generate_quadrature('simpson', 123);
	dm1 = max(diff([0;pts1;1]));
	[pts2, weights] = generate_quadrature('simpson', 121);
	dm2 = max(diff([0;pts2;1]));
	dmax = (dm1+dm2)/2;
	n = quad_points_needed('simpson',dmax);
	%keyboard();
	[status, msg] = assert_eq(n,123);
end
