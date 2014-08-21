%Compute regularisation using Hanson's reg tools library
%
%use a low order approximation to Kd, with trivial quadrature
%for now use the 
function eps = regularise2(Kd, Kdag, Sbar, varargin)
	method = 'lcurve';
	if nargin > 3
		if ischar(varargin{1})
			method = varargin{1};
		end;
	end;
	if nargin > 4
		opts = varargin{2};
	end;
	
	eps = 0;
	switch (method)
		%case 'disc'
		%[x_delta, eps] = discrep(U,sig,V,Sbar,delta,0.69*ones(32,1));
		case 'lcurve'
		eps = lcurve_calculate_eps(Kd, Kdag, Sbar);
		case 'lcurve_lpf'
		eps = lcurve_calculate_eps_lpf(Kd, Kdag, Sbar, opts.Wc);
		otherwise
		error('unknown method for selecting regularisation parameter, %s', method);
	end;
end
