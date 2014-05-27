%Discretise the fredholm operator K, by specifying
%the kernel H(k,z), the number of discrete points n
% and a quadrature method
% assume k,z in [0,1] for now
%sampling points in k domain are same as sampling points in z domain, to
%allow K to be iterated
function Kd = discretise_K(H, n, method)
	[pts, weights] = generate_quadrature(method, n);
	Kd = zeros(n,n);
	row = 1;
	for k = pts
		Kd(row,:) = weights .* arrayfun(@(z) H(k,z), pts);
		row = row+1;
	end;
end
