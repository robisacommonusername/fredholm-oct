function Ax = damp_lsqr(A, epsilon, x, nk, nz, t)
	if strcmp(t,'transp')
		n = length(x)/2;
		Ax = A'*x(1:nk) + epsilon*x((nk+1):(nk+nz));
	else
		Ax = [A*x; epsilon*x];
	end;
end
