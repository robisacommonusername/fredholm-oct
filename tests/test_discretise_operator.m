%Boilerplate
function ret = test_discretise_operator(varargin)
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

function [status, msg] = test_trivial()
	H = @(k,z) exp(-(k^2+z^2));
	[pts,weights] = generate_quadrature('trivial', 2);
	[Kd, Kdag] = discretise_operator(H, pts, weights);
	Ka = 0.5*[exp(-1/8) exp(-5/8); exp(-5/8) exp(-9/8)];
	[status, msg] = assert_eq(Ka, Kd);
end;

function [status, msg] = test_with_A()
	H = @(k,z) exp(-(k^2+z^2));
	[pts,weights] = generate_quadrature('trivial', 2);
	A = [1;2];
	[Kd, Kdag] = discretise_operator(H, pts, weights,A);
	Ka = 0.5*[exp(-1/8) exp(-5/8); 2*exp(-5/8) 2*exp(-9/8)];
	[status, msg] = assert_eq(Ka, Kd);
end;

function [status, msg] = test_gauss_legendre3()
	%use a kernel not symmetric in k, z
	H = @(k,z) k+2*z;
	[pts, weights] = generate_quadrature('gaussn',3);
	[Kd,Kdag] = discretise_operator(H,pts,weights);
	r35 = sqrt(3/5);
	Ka = [1.5*(1-r35), 1+0.5*(1-r35), 0.5*(1-r35)+(1+r35);...
		0.5+(1-r35), 1.5, 0.5+(1+r35);...
		0.5*(1+r35)+(1-r35), 0.5*(1+r35)+1, 1.5*(1+r35)]*diag([5/18;4/9;5/18]);
	[status, msg] = assert_eq(Ka, Kd);
end;

function [status, msg] = test_gauss_legendre3_adjoint()
	%use a kernel not symmetric in k, z
	H = @(k,z) k+2*z;
	[pts, weights] = generate_quadrature('gaussn',3);
	[Kd,Kdag] = discretise_operator(H,pts,weights);
	r35 = sqrt(3/5);
	Ka = transpose([1.5*(1-r35), 1+0.5*(1-r35), 0.5*(1-r35)+(1+r35);...
		0.5+(1-r35), 1.5, 0.5+(1+r35);...
		0.5*(1+r35)+(1-r35), 0.5*(1+r35)+1, 1.5*(1+r35)])*diag([5/18;4/9;5/18]);
	[status, msg] = assert_eq(Ka, Kdag);
end;


