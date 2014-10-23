%Solve the damped least squares problem using Paige and Saunders LSQR
%algorithm, but with a low pass filtering stage

function [chi, error, iterations] = solve_lsqr_lpf(Kd, S, epsilon, filter, varargin)
	%set up options
	if nargin > 4
		opts = varargin{1};
	else
		opts = solve_iteratively_opts(); %defaults
	end;
	x0 = opts.x0;
    if x0 == 0
		[r,c] = size(Kd);
        x0 = zeros(c,1);
    end;
    
    nk = length(S);
    nz = length(x0);
	
	y_damped = [S;sqrt(epsilon)*x0];
	f = @(x,t) filtered_damp_lsqr(Kd, sqrt(epsilon), filter, x, nk, nz, t);
	[chi, flag, error, iterations] = lsqr(f, y_damped, opts.tol, opts.max_iters);
	
	if flag > 0
		warning('LSQR method did not converge');
	end;

end

function Ax = filtered_damp_lsqr(A, epsilon, P, x, nk, nz, t)
	fn_handle = 0;
	if is_function_handle(P)
		f = P;
	else
		f = @(x) P*x;
	end;
	if (strcmp(t,'transp'))
		Ax = f(A'*x(1:nk)) + epsilon*x((nk+1):(nk+nz));
	else
		Ax = [A*f(x); epsilon*x];
	end
end
