%Non adaptive integrator based on Gauss-Legendre 10 point quadrature
%Shouldn't be used directly, designed as a helper function for
%dbl_integrate_adaptive. Unfortuantely due to the vagaries of Matlab/octave
%anonymous and private functions, this has to be in a separate file.
function I = dbl_integrate(f,xa,xb,ya,yb)
	[pts,weights] = generate_quadrature('gauss10',10); %these have been transformed for interval [0,1]
	pts_x = (xb-xa)*pts + xa;
	pts_y = (yb-ya)*pts + ya;
	weights_x = weights*(xb-xa);
	weights_y = weights*(yb-ya);
	[x,y] = meshgrid(pts_x, pts_y);
	I = weights_y'*arrayfun(f,x,y)*weights_x;
end;
