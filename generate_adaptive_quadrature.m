%Generate a set of quadrature points and weights so that a given function
%f can be integrated to a specified tolerance over the interval [0,1].
%Unlike generate_quadrature the interval [0,1] is not subdivided evenly;
% more points are used where the function changes rapidly
%
%params:
%f: function handle of a single variable
%method: string, same as in generate_quadrature
%tol0: absolute tolerance
%pts_per_interval: what is says. Not determined by method, for example
% you can specify method='gauss10' and points_per_interval=20
%max_intervals: maximum number of sub-intervals to divide [0,1] into
function [pts, weights] = generate_adaptive_quadrature(f, method, tol0, pts_per_interval, max_intervals)
	pts_weights = zeros(pts_per_interval*max_intervals, 2);
	pt_offset = 0;
	
	intervals = 1;
	stack = zeros(max_intervals,3);
	sp = 1;
	stack(sp,1) = 0;
	stack(sp,2) = 1;
	stack(sp,3) = tol0;
	
	while sp > 0 && intervals < max_intervals
		a = stack(sp,1);
		b = stack(sp,2);
		tol_n = stack(sp,3);
		sp = sp-1;
		
		[these_pts, these_weights] = generate_quadrature(method, pts_per_interval);
		%transform points and weights from interval [0,1] to [a,b]
		trans_pts = these_pts*(b-a)+a;
		trans_weights = these_weights*(b-a);
		exact = quad(f, a, b, [tol_n/4,0]);
		guess = trans_weights' * arrayfun(f, trans_pts);
		if abs(guess - exact) < tol_n
			pts_weights(pt_offset+1:pt_offset+pts_per_interval,1) = trans_pts;
			pts_weights(pt_offset+1:pt_offset+pts_per_interval,2) = trans_weights;
			pt_offset = pt_offset + pts_per_interval;
		else
			intervals = intervals+2;
			sp = sp+1;
			stack(sp,1) = a;
			stack(sp,2) = (a+b)/2;
			stack(sp,3) = tol_n/2;
			sp = sp+1;
			stack(sp,1) = (a+b)/2;
			stack(sp,2) = b;
			stack(sp,3) = tol_n/2;
		end;
	end;
	
	if sp == 0
		%success
		pts_weights_sorted = sortrows(pts_weights(1:pt_offset,:),1);
		pts = pts_weights_sorted(:,1);
		weights = pts_weights_sorted(:,2);
	else
		warning('integral did not converge, falling back to non-adaptive quadrature');
		[pts,weights] = generate_quadrature(method,pts_per_interval*max_intervals);
	end;
end
