%Calculates the derivative dz/dzeta at a series of points in the zeta
%domain
function dzdzeta = warp_derivative(method,zeta,varargin)
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
			%z = a + (b-a)*zeta
			dzdzeta = (b-a)*ones(length(zeta),1);
		else
			error('Insufficient parameters for linear warping');
		end;
		
	case 'atan'
		if nargin > 2
			alpha = varargin{1};
			%zeta = 2/pi*atan(z/alpha) => z = alpha*tan(pi*zeta/2)
			dzdzeta = alpha*pi/2*(sec(zeta*(pi/2))).^2;
		else
			error('Insufficient parameters for atan warping');
		end;
		
	otherwise
		error('Unrecognised warping method');
	end;
end
