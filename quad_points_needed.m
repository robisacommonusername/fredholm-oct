%Determine how many quad points are needed in interval [0,1] so that no
%two quadrature points are more than d_max apart
%method is a string specifying the quadrature method

function n = quad_points_needed(method,d_max)
	switch(method)
		case 'gauss10'
		max_sep = 0.1488743390;
		%For documentation, here's how to calculate max_sep and what it
		% means:
		%
		%x = [-0.9739065285;-0.8650633667;-0.6794095683;-0.4333953941;...
		%	-0.1488743390;0.1488743390;0.4333953941;0.6794095683;...
		%	0.8650633667;0.9739065285];
		%interval_pts = 0.5*([-2+x(end);x;2+x(1)]+1);
		%max_sep = max(diff(interval_pts));
		%
		%max_sep corresponds to the gap between the 5th and 6th quad points
		num_intervals = ceil(max_sep/d_max);
		n = 10*num_intervals;
		
		case 'simpson'
		n = ceil(1/d_max)+1;
		%ensure odd number of points
		if (mod(n,2) == 0)
			n = n+1;
		end;
		
		case 'trvial'
		n = ceil(1/d_max);
		
		otherwise
		disp('WARNING: Solving by bisection. This may be slow');
		n_min = 888; %initialize, dummy value
		n_max = 1;
		%double n until we bracket the solution
		sep = d_max+1; %dummy init value
		iters = 0;
		max_iters = 100;
		while (sep > d_max && iters < max_iters)
			n_min = n_max;
			n_max = n_min*2;
			[pts,weights] = generate_quadrature(method, n_max);
			sep = max(diff([0;pts;1]));
			iters = iters+1;
		end;
		
		%Now solve by bisection
		dn = n_max-n_min;
		while (dn > 1)
			mid = ceil((n_max+n_min)/2);
			[pts,weights] = generate_quadrature(method, mid);
			sep = max(diff([0;pts;1]));
			if (sep < d_max)
				n_max = mid;
			else
				n_min = mid;
			end;
			dn = n_max-n_min;
		end;
		n1 = min([n_max,n_min]);
		n2 = max([n_max,n_min]);
		[pts_low, weights] = generate_quadrature(method,n1);
		sep = max(diff([0;pts_low;1]));
		if (sep < d_max)
			n = n1;
		else
			n = n2;
		end;
	end;
end
