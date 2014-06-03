%Calculate regularisation parameter, and return Tikhonov regularised
%operator
%Need to specify the quadrature weights so that we can compute |K|
%[Keps, eps] = regularise(K, weights, {method}, {delta}
function [Keps, eps] = regularise(K, varargin)
	method = 'discrepancy';
	delta = 10;
	
	if nargin > 2
		method = varargin{1};
	end;
	if nargin > 3 && method == 'discrepancy'
		delta = varargin{2};
	end;
	
	switch (method)
		case 'discrepancy'
			eps = operator_norm(K,weights)/(1 + delta);
			
		case 'lcurve'
			disp('L curve method not currently implemented');
			eps = 0;
			
		otherwise
			disp('Unkown regularisation method');
			eps = 0;
	end;
	
	[n,n] = size(K);
	Keps = eps*eye(n) + ctranspose(K)*K;
end
