%Inverse of warp data
function z = unwarp_data(method, zeta, varargin)
	switch (method)
	case 'linear'
		if nargin > 2
			a = varargin{1};
			if nargin > 3
				b = varargin{2};
			else
				b = a;
				a = 0;
			end;
			z = a + (b-a)*zeta;
		else
			error('Insufficient parameters for linear warping');
		end;
			
	case 'atan'
		if nargin > 2
			alpha = varargin{1};
			z = alpha*tan(zeta*(pi/2));
		else
			error('Insufficient parameters for atan warping');
		end;
	otherwise
		error('Unrecognised warping method');
	end;
end
