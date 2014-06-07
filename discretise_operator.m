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
	
	%k increases down columns, z increases accross
	%rows, so if x is a column vector, S=Kd*x, S also a col vector
	%hence we need to do [z,k]=meshgrid(...), not [k,z]=mesh...
	[z,k] = meshgrid(pts,pts);
	H_kz = diag(A)*arrayfun(H, k, z);
	W = diag(weights);
	Kd = H_kz*W;
	Kdag = ctranspose(H_kz)*W;
end
