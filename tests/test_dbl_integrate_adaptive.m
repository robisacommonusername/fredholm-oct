%Boilerplate
function ret = test_dbl_integrate_adaptive(varargin)
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

function [status,msg] = test_square()
	%a nice smooth function
	f = @(x,y) exp(-1*(x.^2+y.^2));
	I1 = dbl_integrate_adaptive(f,-1,1,-1,1, 0.0001);
	I2 = dblquad(f,-1,1,-1,1);
	[status,msg] = assert_eq(I1,I2,0.0001,'abs');
end

function [status,msg] = test_rect()
	%a nice smooth function
	f = @(x,y) exp(-1*(x.^2+y.^2));
	I1 = dbl_integrate_adaptive(f,0,3,-1,1, 0.0001);
	I2 = dblquad(f,0,3,-1,1);
	[status,msg] = assert_eq(I1,I2,0.0001,'abs');
end
