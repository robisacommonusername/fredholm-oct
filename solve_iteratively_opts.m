%utility function for parsing input options to pass to solve_iteratively
%usage:
%t = solve_iteratively_opts('key1',val1,'key2',val2,...)
%
% Allowed keys are
% 	x0:
%	norm_k:
%	tol:
%	max_iters:
%	basis:


function opts = solve_iteratively_opts(varargin)
	%set up default values
	opts = struct(...
		'x0',0,...
		'norm_k',0,...
		'max_iters',1000,...
		'tol', 0.0001,...
		'basis','spatial',...,
		'regularise',1);
	
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
