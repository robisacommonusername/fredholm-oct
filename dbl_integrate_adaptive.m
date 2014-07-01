%A function for 2D numerical integration over a rectangular region,
%based on an adaptive Gauss-Legendre 10 point quadrature. Very similar
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
%	f: a function handle, accepting two real arguments. That is
%		f: (RxR) -> C
%	xa, xb: lower and upper limits of integration in x dimension
%	ya, yb: lower and upper limits of integration in y dimension
%varargin
% 	tol: absolute tolerance, default = 4 decimal places. Can pass 0 to use default value
% 	max_iters: maximum recursion depth, default is 100. Cannot be greater
%		than 255, due to limitations of octave interpreter stack
%
%See also: dblquad
%
function [I,status] = dbl_integrate_adaptive(f,xa,xb,ya,yb, varargin)
	tol = 0.0001; %absolute tolerance
	max_iters = 100; %maximum recursion depth
	if nargin > 5
		if varargin{1} > 0
			tol = varargin{1};
		end;
	end;
	if nargin > 6
		if varargin{2} < 255
			max_iters = varargin{2};
		else
			max_iters = 254; %we are already on iteration 1, only do 254 more
		end;
	end;
	
	I=0;
	I1 = dbl_integrate(f,xa,xb,ya,yb);
	if max_iters > 1
		mx = (xa+xb)/2;
		my = (ya+yb)/2;
		I2 = dbl_integrate(f,xa,mx,ya,my)+...
			dbl_integrate(f,mx,xb,ya,my) +...
			dbl_integrate(f,xa,mx,my,yb) +...
			dbl_integrate(f,mx,xb,my,yb);
		%Todo: reuse these values in child calls to prevent double calculation
		if abs(I1-I2) < tol || abs(I2) < tol
			I = I2;
			status = 0;
		else
			[I1,stat1] = dbl_integrate_adaptive(f,xa,mx,ya,my,tol/4,max_iters-1);
			[I2,stat2] = dbl_integrate_adaptive(f,mx,xb,ya,my,tol/4,max_iters-1);
			[I3,stat3] = dbl_integrate_adaptive(f,xa,mx,my,yb,tol/4,max_iters-1);
			[I4,stat4] = dbl_integrate_adaptive(f,mx,xb,my,yb,tol/4,max_iters-1);
			status = max([stat1,stat2,stat3,stat4]);
			I = I1+I2+I3+I4;
		end;
	else
		%recursion limit reached
		I = I1;
		status = 1;
	end;
end;