%Boilerplate
function ret = test_discretise_function(varargin)
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

function [status, msg] = test_undersample_gauss10()
	t=(0:0.01:1)';
	Sexp = cos(10*t);
    [pts,weights] = generate_quadrature('gauss10',20);
	Sd = discretise_function(Sexp, pts, t);
	%pick 50 random test points
	test_points = sort(rand(1,50));
	[status, msg] = assert_eq(interp1(t,Sexp,test_points,'linear','extrap'), interp1(pts,Sd,test_points,'linear','extrap'),0.05);
end

function [status, msg] = test_oversample_gauss10_low_noise()
	t=(0:0.1:1)';
	pure = cos(10*t);
	pwr = 10*log10(trapz(t,pure.^2));
	Sexp = awgn(pure, 25, pwr);
    [pts,weights]=generate_quadrature('gauss10',50);
	Sd = discretise_function(Sexp, pts, t);
	%pick 50 random test points
	test_points = sort(rand(1,50));
	[status, msg] = assert_eq(interp1(t,Sexp,test_points,'linear','extrap'), interp1(pts,Sd,test_points,'linear','extrap'),0.05);
end

function [status, msg] = test_oversample_gauss10_high_noise()
	t=(0:0.1:1)';
	pure = cos(10*t);
	pwr = 10*log10(trapz(t,pure.^2));
	Sexp = awgn(pure, 5, pwr);
    [pts,weights]=generate_quadrature('gauss10',50);
	Sd = discretise_function(Sexp, pts, t);
	%pick 50 random test points
	test_points = sort(rand(1,50));
	[status, msg] = assert_eq(interp1(t,Sexp,test_points,'linear','extrap'), interp1(pts,Sd,test_points,'linear','extrap'),0.05);
end

function [status, msg] = test_undersample_simpson()
	t=(0:0.01:1)';
	Sexp = sin(7*t);
    [pts,weights]=generate_quadrature('simpson',23);
	Sd = discretise_function(Sexp, pts, t);
	%pick 50 random test points
	test_points = sort(rand(1,50));
	[status, msg] = assert_eq(interp1(t,Sexp,test_points,'linear','extrap'), interp1(pts,Sd,test_points,'linear','extrap'),0.05);
end

function [status, msg] = test_oversample_simpson_high_noise()
	t=(0:0.1:1)';
	pure = cos(10*t);
	pwr = 10*log10(trapz(t,pure.^2));
	Sexp = awgn(pure, 5, pwr);
    [pts,weights] = generate_quadrature('simpson',53);
	Sd = discretise_function(Sexp, pts, t);
	%pick 50 random test points
	test_points = sort(rand(1,50));
	[status, msg] = assert_eq(interp1(t,Sexp,test_points,'linear','extrap'), interp1(pts,Sd,test_points,'linear','extrap'),0.05);
end
