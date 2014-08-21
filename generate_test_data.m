%Create dummy interferometric data in 1 dimension
%DEPRECTATE THIS FUNCTION??
%
%params
%Return values
%Sexp: "experimentally" generated interferometric data, as a function of k (fixed Q)
%k_i: sampling points in k domain
%
%inputs
%n: number of points to use
%chi: function handle describing variation of chi sown the z axis
%f: fastcall kernel function
%A: spectrum envelope, sampled in k domain
%A_ki: sampling points in k domain
%alpha: warping parameter
%snr_db: signal to noise ratio in db
%
%use simpson quadrature, as this gives us equally spaced sampling points (more realistic)

function [Sexp, k_i] = generate_test_data(chi, f, fast A, A_ki, alpha, snr_db)
	[pts, weights] = generate_quadrature('simpson',n);
	%use linear remapping
	chibar = chi(alpha*pts);
	k_min = A_ki(1);
	k_max = A_ki(end);
	
	[Hbar, kbar, z] = warp_variables(H, k_min, k_max, alpha);
	Abar = discretise_function(A, pts, kbar(A_ki));
	[Kd,Kdag] = discretise_operator(Hbar, pts, weights, Abar);
	
	Sbar = Kd*chibar;
	k_i = k_min+pts*(k_max-k_min);
	
	%add white noise
	energy = 10*log10((k_max-k_min)* weights' * abs(Sbar).^2);
	Sexp = awgn(Sbar, snr_db, energy);
end
