%utility function for parsing input options to pass to solve_1d
%usage:
%t = solve_1d_opts('key1',val1,'key2',val2,...)
%
% Allowed keys are
%	solver: specify an iterative solver for the second kind equation.
%		Recognised values:
%		'richardson_lpf': default, richardson iteration with projection
%			into subspace of low pass functions
%		'richardson_zero': as above, but without the projection. Vulnerable
%			to noise
%		'richardson_w2': w2 regularisation
%		'bicg_galerkin': galerkin method, using biconjugate gradient solver
% 	n:
%	mean_chi:
%	quad_method:
%	max_iters:
%	tol:
%	basis: specify basis in which kernel will be computed.
%		Recognised values:
%		'spatial': default, basically a Nystrom solver
%		'frequency': use with galerkin method solver when kernel is
%		expressed in basis of complex exponentials
%	correction: correction factors required when solving 1D laminar case
%		(i.e. only solving the Q=0 case). For a Gaussian beam with Ralston's
%		non-uinitary conventions, the required corrections is alpha/k
%	min_feature: minimum feature size. Set to 0 for diffraction limit

function opts = solve_1d_opts(varargin)
	%set up default values
	opts = struct(...
		'solver', 'richardson_lpf',...
		'n',0,...
		'mean_chi',0.8125,...
		'quad_method', 'gauss10',...
		'max_iters',30,...
		'tol', 0.00000001,...
		'basis','spatial',...
		'correction',1,...
		'min_feature',0);
	
	key = '';
	ii = 1;
	while ii <= nargin
		if mod(ii,2) == 1
			%todo: check that key is actually in allowed set
			key = varargin{ii};
		else
			val = varargin{ii};
			opts = setfield(opts, key, val);
		end;
		ii = ii+1;
	end;
end
