%Yet another constrained solver - this solver ensures that the solution
%is restricted to w^2(B). To do this, we do a change of basis on Kd, Kdag
%to express them in a basis of low frequency complex exponentials.
%More correctly, we express chi in terms of the basis of complex exponentials
%but S remains in the k-space domain. i.e K will have different input
%and output bases, i.e. it does an implicit change of basis
%This can be done relatively quickly using an fft

%
%TODO:
% Preconditioning.
%	A possible strategy here is to find the first few eigenfunctions, shrink by appropriate factor, leave other directions unscaled

function [chi, error, iterations] = solve_bicg_galerkin(Kd, S, ki, pts, weights, wc, varargin)
	if nargin > 6
		opts = varargin{1};
	else
		opts = solve_iteratively_opts(); %defaults
	end;
	x0 = opts.x0;
    if x0 == 0
        x0 = zeros(length(pts),1);
    end;
	
	%Choose discretisation level for frequency domain representation of
	%chi based on Nyquist rate
	Ts = pi/wc;
	N = ceil(1/Ts); %recall that length of non-dimensionalised domain is always 1
	%For conveneience, always use an odd length sequence to represent chi
	if (mod(N,2) == 0)
		N = N+1;
	end;
	Ts = 1/N;

	%Convert initial estimate into the frequency domain. Need it to be
	%a sequence of length N, so will need to interpolate
	new_z = transpose(linspace(0,1,N));
	x0_N = interp1(pts, x0, new_z, 'linear', 'extrap');
	x0f = fft(x0_N);
	
	%The frequencies represented are (in order returned by matlab fft)
	%1/(N)*[0, 1/Ts, 2/Ts, ... (N-1)/(2Ts), (1-N)/(2Ts), ... -2/Ts, -1/Ts]
	%Note that N*Ts = 1, since we're working over interval [0,1]
	f = [0:((N-1)/2), ((1-N)/2):-1]';
	[ff,zz] = meshgrid(f,pts);
	exponentials = 1/N*exp(i*2*pi*ff.*zz);
	%If Kd isn't in the right basis, we need to convert it to frequency domain.
	%Note that only chi is in frequency domain, S remains in the k domain
	switch (opts.basis)
		case 'spatial'
		
		%We now need to compute the action of Kd on every exponential of form
		%exp(i*2*pi*f) for the frequencies f defined above
		%This integral will possibly be quite inaccurate due to oscillating integrand. 
		%However, we assume that we have many more quad points than required by Nyquist rate
		%It would be better if the user actually passed a frequency domain
		%kernel directly, since these can be calculated using Filon's
		%rule, etc much more accurately (rather than passing a spatial
		%kernel that need to be converted)
		Kf = Kd*exponentials; 
		
		case 'frequency'
		Kf = Kd;
		
		otherwise
		error('Unrecognised basis/representation for K: "%s"', opts.basis);
	end;
	Kfdag = Kf';
	%Compute regularisation parameter
	if opts.regularise
		epsilon = lcurve_calculate_eps(Kf, Kfdag, S);
	else
		epsilon = 0;
	end;

	%Now solve the system (epsI+Kdagf*Kf)*(chi_tilde) = Kdagf*Sf + eps*x0f
	%Solve the normal equations using conjugate gradient method. Don't
	%Form grammian directly, calculate as Kf'*(Kf*x)
	[chi_tilde,flag,error,iterations] = cgs(@(x) epsilon*x+Kfdag*(Kf*x), Kfdag*S+epsilon*x0f, opts.tol, opts.max_iters);
	if flag > 0
		warning('Conjugate gradient method stagnated or exceeded maximum number of iterations');
	end;
	%Now invert the fourier transform by evaluating the 
	%trig interpolating function at the quadrature points
	chi = exponentials*chi_tilde;
	
	
end
