%A modified version of the lcurve regularisation, that only takes into
%account low frequency noise. i.e, instead of trading off |Ax-b| vs |x|,
%we trade of |Ax-b| vs |Hx|, where H is a low pass filter.
%We will use the solver to force the solution to not have any high freq
%components, rather than try to do this through the regularisation
%Utility function to calculate optimal regularisation parameter using 
%lcurve method. Based on Hansen's code, but doesn't actually plot the
%
% Section copyright Per Christian Hansen
%Wc is the normalised cutoff frequency, it should be < pi, but positive
%We assume that the quadrature points are evenly spaced - this is ok,
%we're using either a trivial or simpson quadrature
function regu = lcurve_calculate_eps_lpf(Kd, Kdag, Sbar, Wc)
	%This is ok, as we're using a low order discretisation
	
	
	%We need to calculate the left and right singular vectors manually
	%rather than just using svd, since Kd* != Kdag in general (i.e when
	%the quad points are not equally spaced). Ok for trivial quadrarure
	%though.
	%[U,s2] = eig(Kd*Kdag);
	%[V,s2] = eig(Kdag*Kd);
	%sm = sqrt(s2);
	%Need to sort the singular values and singular vectors
	[U,sm,V] = csvd(Kd);
	b = Sbar;
	
	npoints = 200;  % Number of points on the L-curve
	smin_ratio = 16*eps;  % Smallest regularization parameter. eps is machine epsilon
	
	[m,n] = size(U);
	[p,ps] = size(sm);
	beta = U'*b; % <U,b> for all the U vectors simultaneously
	beta2 = norm(b)^2 - norm(beta)^2;
	if (ps==1)
		s = sm; beta = beta(1:p);
	else
		s = sm(p:-1:1,1)./sm(p:-1:1,2); beta = beta(p:-1:1);
	end;
	xi = beta(1:p)./s; %unregularised solution
	%Work out cutoffs for filter
	len = length(xi);
	Nc = ceil(Wc/pi * len); %cutoff sample
	upper = ceil(len/2) + Nc;
	lower = ceil(len/2)-Nc;
	if (upper > len) upper = len; end;
	if (lower < 1) lower = 1; end;
		
	eta = zeros(npoints,1); rho = eta; %eta = |Hx|, rho = |Ax-b|
	s2 = s.^2; %singular values squared
	smallest_reg_param = max([s(p),s(1)*smin_ratio]);
	reg_param = transpose(logspace(log10(s(1)), log10(smallest_reg_param), npoints));
	for i=1:npoints
		f = s2./(s2 + reg_param(i)^2);
		spectrum = fftshift(fft(f.*xi));
		eta(i) = norm(spectrum(lower:upper));
		rho(i) = norm((1-f).*beta(1:p)); 
	end;
	if (m > n & beta2 > 0)
		rho = sqrt(rho.^2 + beta2);
	end;
	
	%find corner
	[regu,rho_c,eta_c] = l_corner(rho,eta,reg_param,U,sm,b,'Tikh');
	
end
