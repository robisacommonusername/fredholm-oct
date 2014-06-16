%Similar to discretise K, but here we are resampling our experimental
%data S. Can specify a function for S or a vector. It is assumed that the 
%elements of S are sampled uniformally in the k space, but if not, can 
%specify the current sampling points.
%Just use linear interp. for now.
%assume k in (0,1)

function Sd = discretise_function(Sexp, new_pts, varargin)
	if nargin > 3
		sampling_points = varargin{1};
		if nargin > 4
			%user has specified a warping function kbar(k). Thus we warp
			%the sampling points ki to kbar(ki). Note that the warping
			%function must be able to be called on a vector. This is the
			%case for the function returned from the warp_variables
			%function
			kbar = varargin{2};
			sampling_points = kbar(sampling_points);
		end;
	else
		%assume equally distributed sampling points in [0,1], with
		%half step initial offset if none specified.
		num_points = length(Sexp);
		step = 1.0/num_points;
		sampling_points = step/2:step:(2*num_points-1)*step/2;
	end;
	
	%for prototyping purposes, let the user specify a function
	%to generate the 'experimental' data
	if isa(Sexp, 'function_handle')
		Sexp = arrayfun(Sexp, sampling_points);
	end;
	
	%TODO: low pass filtering to improve noise immunity
	Sd = interp1(sampling_points, Sexp, new_pts, 'linear','extrap');
end
