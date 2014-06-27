%Compute regularisation using Hanson's reg tools library
%
%use a low order approximation to Kd, with trivial quadrature
%for now use the 
function eps = regularise2(H, A, A_ki, S, S_ki, varargin)
	method = 'lcurve';
	if nargin > 5
		if ischar(varargin{1})
			method = varargin{1};
		end;
	end;
	if nargin > 6
		delta = varargin{2};
	elseif strcmp(method,'disc')
		error('no delta parameter specified for use with discrepancy principle');
	end;

	[pts,weights] = generate_quadrature('trivial',32);
	A_resamp = discretise_function(A, pts,A_ki);
	S_resamp = discretise_function(S, pts, S_ki);
	[Kd,Kdag] = discretise_operator(H,pts,weights,A_resamp);
	[U,sig,V] = csvd(Kd);
	
	eps = 0;
	switch (method)
		case 'disc'
		[x_delta, eps] = discrep(U,sig,V,S_resamp,delta,0.69*ones(32,1));
		case 'lcurve'
		eps = lcurve_calculate_eps(Kd, S_resamp);
		otherwise
		error('unknown method for selecting regularisation parameter, %s', method);
	end;
	keyboard();
end
