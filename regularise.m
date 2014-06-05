%Calculate regularisation parameter.

%Need to specify the quadrature weights so that we can compute |K|
%eps = regularise(K, weights, {method}, {delta}
function eps = regularise(K, weights, varargin)
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

end
