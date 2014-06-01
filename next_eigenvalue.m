%Find the next eigenvalue/eigenfunction, given previously known
% eigenfunctions.  Projects discretised operator onto 
%L^2 \setminus \operatorname{span} \{efuncs\}
%efuncs contains the eigenfunctions as a matrix where each eigenfunction
%is a column vector
function [eval, efunc] = next_eigenvalue(Kd, efuncs, varargin)

	%user can optionally specify a tolerance and max number of iterations,
	%otherwise use default values
	maxIters = 1000;
	eps = 0.0001;
	if nargin > 2
		eps = varargin{1};
	end;
	if nargin > 3
		maxIters = varargin{2};
	end;
	
	[rows,cols] = size(Kd);
	x = rand(cols,1);
	x = x/norm(x);
	diff = eps+1;
	guard = 0;
	xold = x;
	efuncst = efuncs';
	
	while (diff > eps) && (guard < maxIters)
		x = Kd*x;
		%project x into the required subspace
		x = x - efuncs*x*efuncst;
		x = x/norm(x);
		diff = norm(x-xold);
		xold = x;
		guard = guard + 1;
	end;
	
	eval = (Kd*x)'*x;
	efunc = x;
end
