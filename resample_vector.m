%This function is designed to replace discretise_function, in  that
%it actually does something useful (rather than just being a wrappper 
%around interp1)
%Use it to compute (for example) A or S at the quadrature points, by
%specifying the quadrature points in the zeta/kappa domain, the warping
%method, etc.  Note that this function is only ever used on k domain data
%which is always warped with the linear transformation
function sigma = resample_vector(S, ki, kappa)
	kmin = ki(1);
	kmax = ki(end);
	if kmin > kmax
		S = flipdim(S)
		ki = flipdim(ki);
		kmin = ki(1);
		kmax = ki(end);
	end;
	
	sigma = interp1((ki - kmin)/(kmax-kmin), S, kappa, 'linear','extrap');
end
