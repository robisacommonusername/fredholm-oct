%Utility function to calculate optimal regularisation parameter using 
%lcurve method. Based on Hansen's code, but doesn't actually plot the
%l-curve. Just copy-paste for now, make readable later
%Always uses tikhonov regularisation
%
% Section copyright Per Christian Hansen
function regu = lcurve_calculate_eps(Kd, Sbar)
	%This is ok, as we're using a low order discretisation
	[U,sm,V] = csvd(Kd);
	b = Sbar;
	
	npoints = 200;  % Number of points on the L-curve
	smin_ratio = 16*eps;  % Smallest regularization parameter. eps is machine epsilon
	
	[m,n] = size(U);
	[p,ps] = size(sm);
	beta = U'*b;
	beta2 = norm(b)^2 - norm(beta)^2;
	if (ps==1)
		s = sm; beta = beta(1:p);
	else
		s = sm(p:-1:1,1)./sm(p:-1:1,2); beta = beta(p:-1:1);
	end;
	xi = beta(1:p)./s;
	
	eta = zeros(npoints,1); rho = eta;
	s2 = s.^2;
	smallest_reg_param = max([s(p),s(1)*smin_ratio]);
	reg_param = transpose(logspace(log10(s(1)), log10(smallest_reg_param), npoints));
	for i=1:npoints
		f = s2./(s2 + reg_param(i)^2);
		eta(i) = norm(f.*xi);
		rho(i) = norm((1-f).*beta(1:p));
	end;
	if (m > n & beta2 > 0)
		rho = sqrt(rho.^2 + beta2);
	end;
	
	%find corner
	[regu,rho_c,eta_c] = l_corner(rho,eta,reg_param,U,sm,b,'Tikh');
	
end