%Discretise the fredholm operator K, by specifying
%the kernel H(k,z), the number of discrete points n
% and a quadrature method
% assume k,z in [0,1] for now
%sampling points in k domain are same as sampling points in z domain, to
%allow K to be iterated

%optional parameter for source spectrum

function [Kd, Kdag] = discretise_operator(H, pts, weights, varargin)
	if nargin > 3
		%A vector or function
		if isa(varargin{1},'function_handle')
			A_func = varargin{1};
			A = A_func(pts);
		else
			A = varargin{1};
		end;
	else
		A = ones(1,length(pts));
	end;

	[k,z] = meshgrid(pts,pts);
	H_kz = diag(A)*arrayfun(H, k, z);
	W = diag(weights);
	Kd = W*H_kz;
	Kdag = W*ctranspose(H_kz);
end
