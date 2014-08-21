%Set quadrature and warping options for fastcall kernels
function opts = fastcall_opts(varargin)
	%set up default values
	opts = struct(...
		'quad_method','gauss10',...
		'n',300,...
		'warp_method','linear',...
		'quad_method_low','trivial',...
		'n_low',64);
	
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
