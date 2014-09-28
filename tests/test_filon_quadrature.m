%Boilerplate
function ret = test_filon_quadrature(varargin)
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

function [status, msg] = test_k0_sin_square()
	%In this case, we reduce down to the standard simpson rule (i.e. no
	%osciallatory component)
	[pts, weights] = filon_quadrature(41, [0]);
	[status, msg] = assert_eq(weights*(sin(2*pi*pts).^2), 0.5);
end

function [status, msg] = test_k0_exp()
	[pts, weights] = filon_quadrature(41,[0]);
	[status, msg] = assert_eq(weights*exp(pts), exp(1)-1);
end

function [status, msg] = test_multik_exp_lowk()
	%Functions of the form exp(i*k*z)*exp(z)
	%Exact solution is 1/(1+ik) * (exp(1+ik) - 1)
	k = (1:20)';
	[pts, weights] = filon_quadrature(81, k);
	exact = 1./(1+i*k).*(exp(1+i*k)-1);
	[zz,kk] = meshgrid(pts, k);
	approx = sum(exp(zz).*weights, 2); %sum across rows (dimension 2)
	[status, msg] = assert_eq(exact, approx);
end


function [status, msg] = test_multik_exp_highk()
	%Functions of the form exp(i*k*z)*exp(z)
	%Exact solution is 1/(1+ik) * (exp(1+ik) - 1)
	k = (10:10:200)';
	[pts, weights] = filon_quadrature(41, k);
	exact = 1./(1+i*k).*(exp(1+i*k)-1);
	[zz,kk] = meshgrid(pts, k);
	approx = sum(exp(zz).*weights, 2); %sum across rows (dimension 2)
	[status, msg] = assert_eq(exact, approx);
end
