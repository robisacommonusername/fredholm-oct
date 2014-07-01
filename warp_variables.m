%Takes a continuous kernel function defined on the rectangle
%[ka,kb]x[0,Inf] and perform a variable transformation
%to get the new kernel
%H'(k', z'), defined on [0,1]x[0,1] with the property
% S(k') = \int_0^1 H'(k',z')\chi(z') dz'
%Also return function corresponding to the variable transformations
%So that S (k') can be evaluated => solve for \chi(z') => get \chi(z)
function [Hbar, kbar, z] = warp_variables(H, ka, kb, zf, varargin)
	method = 'linear';
	gamma = 50;
	if nargin > 4
		method = varargin{1};
	end;
	if nargin > 5
		gamma = varargin{2};
	end; 
	
	%check ka <= kb
	if ka > kb
		temp = ka;
		ka = kb;
		kb = temp;
	end;
	
	%always do linear rescaling on k
	kbar = @(k) (k-ka)/(kb-ka);
	switch(method)
		case 'atan'
		%remap half infinite interval to [0,1] using zbar = 2/pi*atan(z/zf)
		%Will introduce a pole at zbar = 1 due to the dz/dzbar ~ sec^2(pi*zbar/2)
		%This remapping is therefore numerically unstable, and in practice
		%a smooting term should be introduced. Use atan_smooth
		z = @(zbar) zf*tan(zbar*pi/2);
		Hbar = @(ki, zi) H((kb-ka)*ki+ka, zf*tan(zi*pi/2))*zf*pi*0.5.*sec(zi*pi/2).^2;
		
		case 'linear'
		%cut off integration at zf, then rescale linearly
		z =  @(zbar) zbar*zf;
		Hbar = @(ki,zi) H((kb-ka)*ki+ka, zi*zf)*zf;
		
		case 'atan_smooth'
		%TODO: maybe we need to rescale k by a factor of 1/2, we're wasting
		%points past zbar = 0.5
		
		%same as atan remapping, but introduce an exponential smoothing
		%term to force chi to zero as z->infinity. Smooths out the poles
		%introduced by dz/dzbar ~ sec^2(zbar)
		z = @(zbar) zf*tan(zbar*pi/2);
		Hbar = @(ki, zi) H((kb-ka)*ki+ka, z(zi))/(1+exp(gamma*(z(zi)-zf)))*zf*pi*0.5.*sec(zi*pi/2).^2;
		
		otherwise
		disp('WARNING: unrecognised variable warping method. Falling back to linear');
		z =  @(zbar) zbar*zf;
		Hbar = @(ki,zi) H((kb-ka)*ki+ka, zi*zf)*zf;
		
	end;
	
	
end
