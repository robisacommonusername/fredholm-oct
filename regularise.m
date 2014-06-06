%Calculate regularisation parameter.

%Need to specify the quadrature weights so that we can compute |K|
%args
%method: either 'disc' or 'lcurve'
%if method == disc, params are
%normK = |K|
%delta = SNR

%if method == 'lcurve', args are
%K, Kdag, weights
function eps = regularise(varargin)
	if nargin > 0
		ii = 1;
		if ischar(varargin{1})
			method = varargin{1};
			ii = 2;
		else
			method = 'disc';
		end;
		
		switch (method)
			case 'disc'
			if nargin >= 3
				normK = varargin{ii};
				delta = varargin{ii+1};
				eps = normK/(1+delta);
			else
				disp('WARNING: insufficient parameters for discrepancy prnciple')
				eps = 0;
			end;
			
			case 'lcurve'
			disp('WARNING: lcurve method is not currently implemented');
			eps = 0;
			
			otherwise
			disp('WARNING: unrecognised regularisation method');
			eps = 0;
		end;
	else
		disp('WARNING: no parameters entered for regularisation');
		eps = 0;
	end;

end
