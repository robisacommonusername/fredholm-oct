%Compute a filon quadrature for integrating functions of the form
%exp(i*k*z)f(z) for a range of different (large) k values
%The weights will be different for each different k, so it is a matrix.
%we will integrate over the interval [0,1]
%See Davis and Rabinowitz, "Methods of Numerical Integration", pp. 151-5, 1984
function [pts, weights] = filon_quadrature(npoints, k)
	[r,c] = size(k);
	if r < c
		warning('k should be specified as a column vector');
		k = k';
		r = c;
	end;
	
	if mod(npoints,2) == 0
		warning('Number of evaluation points for Filon quadrature should be odd');
		npoints = npoints + 1;
	end;
	
	%N is the number of "panels", it should be even
	N = npoints-1;
	
	h = 1/N;
	theta = k*h;
	%Compute the alpha, beta, gamma parameters. Use the taylor series for
	%small theta
	smalls = theta < 1/6;
	bigs = theta >= 1/6;
	
	%Compute small parameters
	ts = theta(smalls);
	ts2 = ts.^2;
	ts3 = ts2.*ts;
	alpha_s = ts3.*(2/45 + ts2.*(-2/315 + ts2.*2/4725));
	beta_s = 2/3 + ts2.*(2/15 + ts2.*(-4/105 + ts2.*(2/567)));
	gamma_s = 4/3 + ts2.*(-2/15 + ts2.*(1/210 + ts2.*-1/11340));
	
	%Compute the large parameters
	tb = theta(bigs);
	t2 = tb.^2;
	t3 = t2.*tb;
	sint = sin(tb);
	cost = cos(tb);
	sintcost = sint.*cost;
	alpha_b = (t2 + tb.*sintcost - 2*sint.^2)./t3;
	beta_b = 2*(tb.*(1+cost.^2) - 2*sintcost)./t3;
	gamma_b = 4*(sint-tb.*cost)./t3;
	
	alpha = zeros(r,1);
	alpha(smalls) = alpha_s;
	alpha(bigs) = alpha_b;
	beta = zeros(r,1);
	beta(smalls) = beta_s;
	beta(bigs) = beta_b;
	gamma = zeros(r,1);
	gamma(smalls) = gamma_s;
	gamma(bigs) = gamma_b;
	
	
	pts = transpose(linspace(0,1,npoints));
	pts_even = pts(1:2:npoints);
	pts_odd = pts(2:2:N);
	%let k go down column, weights will extend across rows
	[zz_e,kk_e] = meshgrid(pts_even, k);
	[zz_o,kk_o] = meshgrid(pts_odd, k);
	
	C_e = exp(i*kk_e.*zz_e);
	C_e(:,1) = 0.5*C_e(:,1);
	C_e(:,end) = 0.5*C_e(:,end);
	C_o = exp(i*kk_o.*zz_o);
	
	weights = zeros(r, npoints);
	weights(:,1:2:npoints) = diag(beta)*C_e;
	weights(:,2:2:N) = diag(gamma)*C_o;
	
	%Add in the extra terms at z=0 and z=1
	weights(:,1) = weights(:,1) + i*alpha;
	weights(:,end) = weights(:,end) - i*alpha.*exp(i*k);
	
	%Now scale by h
	weights = h*weights;
	
end
