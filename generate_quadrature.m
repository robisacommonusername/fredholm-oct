%generate the evaluation points and weights
%for various quadrature methods. Assume
%integral is over [0,1]
function [pts, weights] = generate_quadrature(method, n)
	switch(method)
		case 'gauss10'
		%lookup
		w = [0.0666713443;0.1494513492;0.2190863625;0.2692667193;...
			0.2955242247;0.2955242247;0.2692667193;0.2190863625;...
			0.1494513492;0.0666713443];
			
		x = [-0.9739065285;-0.8650633667;-0.6794095683;-0.4333953941;...
			-0.1488743390;0.1488743390;0.4333953941;0.6794095683;...
			0.8650633667;0.9739065285];
			
		%classic gauss-legendre 10 point quadrature. Subdivide interval as necessary
		num_intervals = ceil(n/10);
		if (mod(n,10) ~= 0)
			warning('order (n) should be a multiple of 10 when using Gauss 10 point quadrature');
		end;
		
		interval_width = 1.0/num_intervals;
		weights = repmat(interval_width/2*w,num_intervals,1);
		pts = zeros(n,1);
		for interval_i = 1:num_intervals
			a = (interval_i-1)*interval_width;
			b = interval_i*interval_width;
			pts(10*(interval_i-1)+1:10*(interval_i-1)+10) = 0.5*(b+a)+0.5*(b-a)*x;
		end;
			
		case 'gaussn'
		%n point gauss-legendre quadrature. Evaluation points are the roots
		% of the Legendre polynomials
		
		%1. Calculate the nth Legendre polynomial recursively (Bonnet formula)
		% n P[n](x) = (2*n-1)xP[n-1](x) - (n-1)P[n-2](x)
		Plast = [1];
		Pthis = [1 0];
		for ii = 2:n
			temp = Pthis;
			Pthis = 1/ii * ((2*ii - 1)*conv(Pthis, [1,0]) - (ii-1)*[0 0 Plast]);
			Plast = temp;
		end;
		
		%2. Evaluation points. Remap for interval [0,1]
		xi = sort(roots(Pthis));
		pts = (xi+1)/2;
		
		%3. Weights. The weights are 2/[(1-xi^2)P'(xi)^2] * (1-0)/2
		% factor of 1/2 is due to shifting integration domain to [0,1] from [-1,1]
		Pderiv = polyder(Pthis);
		weights = 1 ./ ((1 - xi.^2) .* (arrayfun(@(x) polyval(Pderiv,x), xi)).^2);
		
		case 'clenshawcurtis10'
		
		case 'clenshawcurtisn'
		%UNTESTED - DO NOT USE
		%n point clenshaw curtis quadrature. Evaluation points
		%are the roots of the first kind Chebyshev polynomials
		%Tn(x) = cos (n acos(x)) = 0
		N = n-1;
		
		
		k = (1:N/2)';
		pts = ([-1*cos(k*pi/N);1;cos(k*pi/N)]+1)/2;
		d = [1;2./(1-(2:2:(N-2))'.^2);1/(1-N^2)];
		w = dct(d)/2;
		weights = [flipud(w);w(2:end)];
		
		case 'simpson'
		%subdivide interval [0,1]
		%because of shared end points, n = 2*num_intervals+1
		num_intervals = ceil((n-1)/2);
		if mod(n,2) ~= 1
			warning('n should be odd for use with Simpsons rule');
		end;
		
		interval_width = 1.0/num_intervals;
		step = interval_width/2;
		pts = (0:step:1)';
		weights = [interval_width/6;repmat([2*interval_width/3;interval_width/3],num_intervals,1)];
		%fix last entry, don't double count
		weights(end) = interval_width/6;
		
		case 'trivial'
		step = 1/n;
		weights = step * ones(n,1);
		pts = (step/2:step:1-step/2)';
		
		case 'filon'
		%Generates the filon weights/abscissas assuming the oscillating
		%term is actually a constant (i.e. k=0). Should be compatible
		%(more or less) with simpson. One difference is that 'simpson'
		%doesn't use end points (i.e. open quadrature), but 'filon' will 
		%return the classic simpson formulae, including the end points
		%(closed quadrature)
		%Do NOT use these wights/pts for
		%discretising the kernel, use them for integrating the slowly
		%varying part (i.e. chi). For discretising the kernel with filon
		%method, use filon_quadrature, and specify the ki
		[pts, weights_transpose] = filon_quadrature(n, 0);
		weights = weights_transpose';
		
		otherwise
		warning('unknown quadrature method. Falling back to trivial method');
		[pts, weights] = generate_quadrature('trivial', n);
		
		
	end;
end
