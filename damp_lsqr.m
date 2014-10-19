function Ax = damp_lsqr(A, epsilon, x, t)
	if strcmp(t,'transp')
		n = length(x)/2;
		Ax = A'*x(1:n) + epsilon*x((n+1):(2*n));
	else
		Ax = [A*x; epsilon*x];
	end;
end
