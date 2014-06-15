%Kernel function for paraxial gaussian beam.
%Returns an annonymous function for given Q, alpha, z0
%Q is just a scalar here, as kernel only depends on $\|mathbf{Q}\|$
function H = gauss_kernel(Q,alpha,z0)
	H = @(k,z) 1/4/pi/(alpha^2/k^2 + i*(z-z0)/k)*exp(2*i*(z-z0)*k-Q^2/2*(alpha^2/2/k^2+i*(z-z0)/2/k));

end
