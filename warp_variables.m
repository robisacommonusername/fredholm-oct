%Takes a continuous kernel function defined on the rectangle
%[ka,kb]x[0,Inf] and perform a variable transformation
%to get the new kernel
%H'(k', z'), defined on [0,1]x[0,1] with the property
% S(k') = \int_0^1 H'(k',z')\chi(z') dz'
%Also return function corresponding to the variable transformations
%So that S (k') can be evaluated => solve for \chi(z') => get \chi(z)
function [Hbar, kbar, z] = warp_variables(H, ka, kb, varargin)
	alpha = 1;
	if nargin > 3
		alpha = varargin{1};
	end;
	
	%check ka <= kb
	if ka > kb
		temp = ka;
		ka = kb;
		kb = temp;
	end;
	
	kbar = @(k) (k-ka)/(kb-ka);
	z = @(zbar) alpha*tan(zbar*pi/2);
	
	Hbar = @(ki, zi) H((kb-ka)*ki+ka, alpha*tan(zi*pi/2))*alpha*pi*0.5*sec(zi*pi/2)^2;
end
