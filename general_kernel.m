%Specify fourier transform of beam profile, and compute value of kernel
%for a given k,z.  An expensive function
%Q is a vector
%No variable remapping
function H = general_kernel(g0,Q,z0)
	%todo - estimate integration interval better
	%axial wavelength helper function - no proper 'let', will be called twice
	aw = @(q,k) sqrt(k.^2-(q*q'));
	H = @(k,z) dbl_integrate_adaptive(@(q1,q2) 1/aw([q1,q2], k)*exp(i*aw([q1,q2],k)*(z-z0))*g0([q1,q2], k)*exp(i*aw(Q-[q1,q2], k)*(z-z0))*g0(Q-[q1,q2], k),-0.5*k,k/2,-0.5*k,k/2,0,5);
end

