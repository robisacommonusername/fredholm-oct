%utility function for parsing input options to pass to solve_1d
%usage:
%t = solve_1d_opts('key1',val1,'key2',val2,...)
%
% Allowed keys are
% 	n:
%	mean_chi:
%	quad_method:
%	max_iters:
%	tol:
%	reg_method:
%	snr_db:
%

function opts = solve_1d_opts(varargin)
	%set up default values
	opts = struct(...
		'n',0,...
		'mean_chi',0.8125,...
		'quad_method', 'gauss10',...
		'max_iters',10000,...
		'tol', 0.0001,...
		'reg_method','lcurve',...
		'snr_db', 10);
	
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
