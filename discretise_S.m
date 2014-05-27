%Similar to discretise K, but here we are resampling our experimental
%data S. Can specify a function for S or a vector. It is assumed that the 
%elements of S are sampled uniformally in the k space, but if not, can 
%specify the current sampling points.
%Just use linear interp. for now.
%assume k in (0,1)
function Sd = discretise_S(Sexp, n, method, varargin)
	if nargin > 3
		sampling_points = varargin{1};
	else
		num_points = rows(Sexp);
		step = 1.0/num_points;
		sampling_points = step/2:step:(2*num_points-1)*step/2;
	end;
	
	%for prototyping purposes, let the user specify a function
	%to generate the 'experimental' data
	if isa(Sexp, 'function_handle')
		Sexp = arrayfun(Sexp, sampling_points);
	end;
	
	new_pts = generate_quadrature(method, n);
	Sd = zeros(n,1);
	
	%resample. Use linear interpolation
	derivs = diff(Sdiff) ./ diff(sampling_points)';
	last = derivs(rows(derivs));
	dS = [derivs; last] %add the last element on twice
	sampling_point_idx = 1;
	num_sampling_points = columns(sampling_points);
	for new_pt_idx = 1:n
		%find the closest current sampling point to the new sampling point
		new_x = new_pts(new_pt_idx);
		closest = sampling_points(sampling_point_idx);
		old_closest = closest;
		this_dist = abs(new_x - closest);
		prev_dist = this_dist + 1;
		
		while (prev_dist > this_dist)
			sampling_point_idx = sampling_point_idx + 1;
			old_closest = closest;
			%need to be careful when we overshoot last sampling point
			if (sampling_point_idx <= num_sampling_points)
				closest = sampling_points(sampling_point_idx);
				prev_dist = this_dist;
				this_dist = abs(new_x - closest);
			end;
		end;
		
		%we will have overshot the closest point by 1
		closest = old_closest;
		sampling_point_idx = sampling_point_idx - 1;
		
		%now perform linear interpolation
		Sd(new_pt_idx) = Sexp(sampling_point_idx) + dS(sampling_point_idx)*(new_x - closest);
	end;
	
end
