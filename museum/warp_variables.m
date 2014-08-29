%Takes a continuous kernel function defined on the rectangle
%[ka,kb]x[0,Inf] and perform a variable transformation
%to get the new kernel
%H'(k', z'), defined on [0,1]x[0,1] with the property
% S(k') = \int_0^1 H'(k',z')\chi(z') dz'
%Also return function corresponding to the variable transformations
%So that S (k') can be evaluated => solve for \chi(z') => get \chi(z)
function [kbar_of_k, k_of_kbar, zbar_of_z, z_of_zbar, dz_dzbar] = ...
	warp_variables(ka, kb, zf, varargin)
	
	method = 'linear';
	gamma = 50;
	if nargin > 3
		method = varargin{1};
	end;
	if nargin > 4
		gamma = varargin{2};
	end; 
	
	%check ka <= kb
	if ka > kb
		temp = ka;
		ka = kb;
		kb = temp;
	end;
	
	%always do linear rescaling on k
	kbar_of_k = @(k) (k-ka)/(kb-ka);
	k_of_kbar = @(kbar) (kb-ka)*kbar + ka;
	switch(method)
		case 'atan'
		%remap half infinite interval to [0,1] using zbar = 2/pi*atan(z/zf)
		%Will introduce a pole at zbar = 1 due to the dz/dzbar ~ sec^2(pi*zbar/2)
		%This remapping is therefore numerically unstable, and in practice
		%a smooting term should be introduced. Use atan_smooth
		z_of_zbar = @(zbar) zf*tan(zbar*pi/2);
		zbar_of_z = @(z) 2/pi*atan(z/zf);
		dz_dzbar = @(zbar) 0.5*zf*pi*(sec(zbar*pi/2)).^2;
		
		case 'linear'
		%cut off integration at zf, then rescale linearly
		z_of_zbar =  @(zbar) zbar*zf;
		zbar_of_z = @(z) z/zf;
		dz_dzbar = @(zbar) zf;
		
		otherwise
		warning('unrecognised variable warping method. Falling back to linear');
		z_of_zbar =  @(zbar) zbar*zf;
		zbar_of_z = @(z) z/zf;
		dz_dzbar = @(zbar) zf;
		
	end;
	
	
end
