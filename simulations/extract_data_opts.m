%utility function for parsing input options to pass to solve_iteratively
%usage:
%t = solve_iteratively_opts('key1',val1,'key2',val2,...)
%
% Allowed keys are
% 	EH: 0 = Extract E field, 1 = extract H field
%	component: component of field to extract: 1=x,2=y,3=z
%	plane: Extract data from the plane (x|y|z) = a. Use x=1,y=2,z=3
%	a: constant for plane on which to extract data
%	out: %Filename to save extracted data (if desired)

function opts = extract_data_opts(varargin)
	%set up default values
	opts = struct(...
		'EH',0,...
		'component',1,...
		'plane',3,...
		'a', 0,
		'out','');
	
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

