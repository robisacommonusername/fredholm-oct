%This function is a lot like warp_variables, except instead of returning
%a function handle, it works on actual data

function zeta = warp_data(method,z,varargin)
	switch (method)
	case 'linear'
			%zeta = (z-a)/(b-a)
			if nargin > 2
				a = varargin{1};
				if nargin > 3
					b = varargin{2};
				else
					b = a;
					a = 0;
				end;
				zeta = (z-a)/(b-a);
			else
				error('Insufficient parameters for linear warping');
			end;
			
	case 'atan'
			if nargin > 2
				alpha = varargin{1};
				zeta = 2/pi*atan(z/alpha);
			else
				error('Insufficient parameters for atan warping');
			end;
			
	otherwise
		error('Unrecognised warping format');
	end;
end
