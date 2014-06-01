%Specify fourier transform of beam profile, and compute value of kernel
%for a given k,z.  An expensive function
%Q is a vector
%No variable remapping
function H = general_kernel(g0,Q,z0)
	%todo - estimate integration interval better
	%axial wavelength helper function - no proper 'let', will be called twice
	aw = @(k,q1,q2) sqrt(k^2-q1^2-q2^2);
	awr = @(k,r) sqrt(k^2-r^2);
	%H = @(k,z) dblquad(@(q1,q2) 1/aw(k,q1,q2)*exp(i*aw(k,q1,q2)*(z-z0))*g0(k,q1,q2)*exp(i*aw(k,Q(1)-q1, Q(2)-q2)*(z-z0))*g0(k,Q(1)-q1,Q(2)-q2),-1*k,k,-1*k,k);
	
	%change to radial coordinates, and use circular region of integration
	H = @(k,z) dblquad(@(r,t) r/awr(k,r)*exp(i*awr(k,r)*(z-z0))*g0(k,r*cos(t),r*sin(t))*exp(i*aw(k,Q(1)-r*cos(t),Q(2)-r*sin(t))*(z-z0))*g0(k,Q(1)-r*cos(t),Q(2)-r*sin(t)),0,k,0,2*pi);
end
