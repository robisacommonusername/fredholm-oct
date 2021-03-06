%As for dbl_integrate_adaptive, but instead of passing a function handle, we use
%an irregular grid of points and function values specified on that grid, and interpolate
%
%A function for 2D numerical integration over a rectangular region,
%based on an adaptive Gauss-Kronrod (G7-K15) point quadrature. Very similar
%to built in function dblquad, but without the restrictions: can be
%be called recursively, can integrate complex functions, and can be
% called on non-vectorised functions.
%
%Usage:
%	[I,status] = dbl_integrate_adaptive(f,xa,xb,ya,yb, varargin)
%
%Outputs:
%	I: the computed integral
%	status: a status code
%		0: success
%		1: recursion limit exceeded
%
%Inputs:
%	xx: matrix of x values on which f is defined (meshgrid format)
%	yy: matrix of y values on which f is defines (meshgrid format)
%	f: 	matrix of function values f(xx,yy)
%	xa, xb: lower and upper limits of integration in x dimension
%	ya, yb: lower and upper limits of integration in y dimension
%varargin
% 	tol: absolute tolerance, default = 4 decimal places. Can pass 0 to use default value
% 	max_iters: maximum recursion depth, default is 10. Cannot be greater
%		than 255, due to limitations of octave interpreter stack. Choosing
%		more than 255 only proves that you are an idiot anyway.
%
%See also: dblquad
%
function [I,status] = dbl_integrate_adaptive_interpolate_irregular(x,y,f,xa,xb,ya,yb, varargin)
	tol = 0.0001; %absolute tolerance
	max_iters = 10; %maximum recursion depth
	if nargin > 7
		if varargin{1} > 0
			tol = varargin{1};
		end;
	end;
	if nargin > 8
		if varargin{2} < 255
			max_iters = varargin{2};
		else
			max_iters = 254; %we are already on iteration 1, only do 254 more
		end;
	end;
	
	%Select x and y that are within range
	indices = x>xa & x<xb & y>ya & y<yb;
	x = x(indices);
	y = y(indices);
	f = f(indices);
	
	[I,err] = dbl_integrate_interpolate_irregular(x,y,f,xa,xb,ya,yb);
	status = 0;
	if (max_iters > 1)
		if (err > tol)
			mx = (xa+xb)/2;
			my = (ya+yb)/2;
			[I1,stat1] = dbl_integrate_adaptive(x,y,f,xa,mx,ya,my,tol/4,max_iters-1);
			[I2,stat2] = dbl_integrate_adaptive(x,y,f,mx,xb,ya,my,tol/4,max_iters-1);
			[I3,stat3] = dbl_integrate_adaptive(x,y,f,xa,mx,my,yb,tol/4,max_iters-1);
			[I4,stat4] = dbl_integrate_adaptive(x,y,f,mx,xb,my,yb,tol/4,max_iters-1);
			status = max([stat1,stat2,stat3,stat4]);
			I = I1+I2+I3+I4;
		end;
	else
		%recursion limit reached
		status = 1;
	end;
end
