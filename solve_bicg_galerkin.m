%Yet another constrained solver - this solver ensures that the solution
%is restricted to w^2(B). To do this, we do a change of basis on Kd, Kdag
%to express them in a basis of low frequency complex exponentials.
%This can be done relatively quickly using an fft
%This method has a rather large setup time in some cases. It is most
%effective when the quadrature points are much more closely spaced than
%the Nyquist sampling rate
%
%TODO:
% Preconditioning.
%	A possible strategy here is to find the first few eigenfunctions, shrink by appropriate factor, leave other directions unscaled

function [chi, error, iterations] = solve_bicg_galerkin(Kd, S, pts, weights, wc, varargin)
	if nargin > 5
		opts = varargin{1};
	else
		opts = solve_iteratively_opts(); %defaults
	end;
	x0 = opts.x0;
    if x0 == 0
        x0 = zeros(length(S),1);
    end;
	
	%Choose discretisation level for frequency domain representation
	%based on Nyquist rate
	Ts = pi/wc;
	N = ceil(1/Ts); %recall that length of non-dimensionalised domain is always 1
	%For conveneience, always use an odd length sequence to represent chi
	if (mod(N,2) == 0)
		N = N+1;
	end;
	
	%Our data is over-sampled - we only need S at N wavenumbers
	k = linspace(0,1,N);
	Sint = interp1(pts,S,k,'spline','extrap'); %We want a fairly smooth interpolation, I think cubic splines are ok
	Sf = transpose(fft(Sint)); %don't use ', that does ctranspose
	x0int = interp1(pts,x0,k,'spline','extrap');
	x0f = transpose(fft(x0int));
	
	%The frequencies represented are (in order returned by matlab fft)
	%1/N*[0, 1/Ts, 2/Ts, ... (N-1)/(2Ts), (1-N)/(2Ts), ... -2/Ts, -1/Ts]
	f = 1/(N*Ts)*[0:((N-1)/2), ((1-N)/2):-1]';
	
	%If Kd isn't in the right basis, we need to convert it to frequency domain
	switch (opts.basis)
		case 'spatial'
		%keyboard();
		%we don't need Kd defined at so many K values - lets cut down
		n_quad = length(pts);
		Kdint = zeros(N,n_quad);
		for ii = 1:n_quad
			Kdint(:,ii) = interp1(pts,Kd(:,ii),k,'spline','extrap');
		end;
		%keyboard();
		%We now need to compute the action of Kd on every exponential of form
		%exp(i*2*pi*f) for the frequencies f defined above
		[ff,zz] = meshgrid(f,pts); %f increases across row, z down columns
		KdE = 1/N*Kdint*exp(i*2*pi*ff.*zz); 
		
		%Now we fourier transform every column of KdE to express it in freq
		%domain. Calling fft on matrix automatically ffts each column
		Kf = fft(KdE);
		
		case 'frequency'
		Kf = Kd;
		
		otherwise
		error('Unrecognised basis/representation for K: "%s"', opts.basis);
	end;
	
	%Compute regularisation parameter
	if opts.regularise
		epsilon = lcurve_calculate_eps(Kf, Kf', Sf);
	else
		epsilon = 0;
	end;

	%Now solve the system (epsI+Kdagf*Kf)*(chi_tilde) = Kdagf*Sf + eps*x0f
	chi_tilde_0 = fft(x0);
	[chi_tilde,flag,error,iterations] = bicg(epsilon*eye(N)+Kf'*Kf, Kf'*Sf+epsilon*x0f, opts.tol, opts.max_iters);
	if flag > 0
		warning('Biconjugate gradient method stagnated or exceeded maximum number of iterations');
	end;
	%Now invert the fourier transform by evaluating the 
	%trig interpolating function at the quadrature points
	chi = 1/N*exp(i*2*pi*ff.*zz)*chi_tilde;
	
	
end
