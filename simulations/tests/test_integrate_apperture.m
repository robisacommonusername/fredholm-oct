function ret = test_integrate_apperture(varargin)
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

function [status, msg] = test_integrate_gaussian()
	%integrate the function exp(-(x^2+y^2)/2) over disc r < 5. Should get approx 2pi
	%We will use a polar grid, which will tend to concentrate the 
	%evaluation points near r=0 (which is what we want anyway)
	r = linspace(0,5,100);
	theta = linspace(0,2*pi,100);
	[rr,tt] = meshgrid(r,theta);
	x = rr.*cos(tt);
	y = rr.*sin(tt);
	field = vec(exp(-0.5*rr.^2));
	vertices = [vec(x),vec(y),zeros(100^2,1)];
	result = integrate_apperture(vertices, field);
	[status,msg] = assert_eq(result,2*pi,0.01);
end

function [status, msg] = test_integrate_many_gaussians()
	%Similar to previous test, but this time we integrate many gaussians
	%changing the "standard deviation" each time
	%integrate the function exp(-(x^2+y^2)/2sigma^2) over disc r < 5. 
	%Should get approx 2*pi*sigma^2
	%We will use a polar grid, which will tend to concentrate the 
	%evaluation points near r=0 (which is what we want anyway)
	r = linspace(0,5,100);
	theta = linspace(0,2*pi,100);
	sigma = linspace(0.5,1.5,5)';
	field = zeros(100^2,5);
	[rr,tt] = meshgrid(r,theta);
	x = rr.*cos(tt);
	y = rr.*sin(tt);
	vertices = [vec(x),vec(y),zeros(100^2,1)];
	for is = 1:5
		field(:,is) = vec(exp(-0.5/sigma(is)^2*rr.^2));
	end;
	
	result = integrate_apperture(vertices, field);
	[status,msg] = assert_eq(result,2*pi*sigma.^2,0.01);
end

function [status, msg] = test_mixed_signs()
	%Use a function that goes above and below the z axis
	%integrate the function exp(-(x^2+y^2)/2)-0.1 over disc r < 5. 
	%Should get approx 2pi - 2.5pi = -0.5*pi
	%We will use a polar grid, which will tend to concentrate the 
	%evaluation points near r=0 (which is what we want anyway)
	r = linspace(0,5,100);
	theta = linspace(0,2*pi,100);
	[rr,tt] = meshgrid(r,theta);
	x = rr.*cos(tt);
	y = rr.*sin(tt);
	field = vec(exp(-0.5*rr.^2)-0.1);
	vertices = [vec(x),vec(y),zeros(100^2,1)];
	result = integrate_apperture(vertices, field);
	[status,msg] = assert_eq(result,-0.5*pi,0.01);
end

function [status, msg] = test_complex()
	%integrate a complex value function
	%integrate the function exp(-(x^2+y^2)/2)+2*i*exp(-(x^2+y^2)/3)
	% over disc r < 5. Should get approx 2pi + i*6pi
	%We will use a polar grid, which will tend to concentrate the 
	%evaluation points near r=0 (which is what we want anyway)
	r = linspace(0,5,100);
	theta = linspace(0,2*pi,100);
	[rr,tt] = meshgrid(r,theta);
	x = rr.*cos(tt);
	y = rr.*sin(tt);
	field = vec(exp(-0.5*rr.^2) + 2*i*exp(-1/3*rr.^2));
	vertices = [vec(x),vec(y),zeros(100^2,1)];
	result = integrate_apperture(vertices, field);
	[status,msg] = assert_eq(result,2*pi+6*i*pi,0.01);
	
end
