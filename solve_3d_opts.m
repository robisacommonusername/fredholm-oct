%utility function for parsing input options to pass to solve_3d
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
%		'bicg_galerkin': galerkin method, using biconjugate gradient solver
%		'cgls': conjugate gradient on normal equations, solution in L^2
%		'cgls_lpf': conjugate gradient on normal equations, solution in w^2(B)
%		'lsqr': paige and saunders lsqr algorithm, L^2 solution
%		'lsqr_lpf': paige and saunders algorithm, solution in w^2(B)
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
%	downsample: if set to 1, will interpolate user supplied data S to
%		only use values of S at the collocation points. By default, this is
%		set to 0, in which case all the data is used (possibly overdetermining
%		the poblem), and a least squares solution returned. The only reason
%		one might wish to downsample is to improve performance, as it will
%		downsample S to the minimum number of entries required to solve the problem
%	do_filter: Contrain solution to w^2(B)? 1=yes, 0=no, default=1

function opts = solve_3d_opts(varargin)
	%set up default values
	opts = struct(...
		'solver', 'lsqr_lpf',...
		'n',0,...
		'mean_chi',0.8125,...
		'quad_method', 'gauss10',...
		'max_iters',0,...
		'tol', 0.000001,...
		'basis','spatial',...
		'correction',1,...
		'min_feature',0,...
		'downsample',0,
		'do_filter',1);
	
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
