%Similar to discretise K, but here we are resampling our experimental
%data S. Can specify a function for S or a vector. It is assumed that the 
%elements of S are sampled uniformally in the k space, but if not, can 
%specify the current sampling points.
%Just use linear interp. for now.
%assume k in (0,1)
function Sd = discretise_function(Sexp, n, method, varargin)
	if nargin > 3
		sampling_points = varargin{1};
	else
		num_points = length(Sexp);
		step = 1.0/num_points;
		sampling_points = step/2:step:(2*num_points-1)*step/2;
	end;
	
	%for prototyping purposes, let the user specify a function
	%to generate the 'experimental' data
	if isa(Sexp, 'function_handle')
		Sexp = arrayfun(Sexp, sampling_points);
	end;
	
	new_pts = generate_quadrature(method, n);
	
	%TODO: low pass filtering to improve noise immunity
	Sd = interp1(sampling_points, Sexp, new_pts, 'linear','extrap');
end
