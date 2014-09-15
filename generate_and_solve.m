%Generate some test data, then attempt to solve the problem. This function
%is largely for testing purposes. This is a 1d function

%Input parameters
%	chi: function handle, describing chi down the z axis
%	f: fastcall function
%	A: spectral envelope
%	ki: sampling points for spectral envelope
%	zf: sample thickness
%	
%	varargin:
%		noise_ratio: noise_power/signal power. Defaults to 0
%		opts: a solve_1d_opts structure
function [chi_exp,z_exp] = generate_and_solve(chi,f,A,ki,zf,varargin)
	%We specify the noise ratio, rather than snr, as this allows us
	%to specify zero noise. The snr in db is -10 log10(noise_ratio)
	noise_ratio = 0;
	if nargin > 5
		noise_ratio = varargin{1};
	end;
	
	%read opts for the 1d solver
	if nargin > 6
		opts = varargin{2};
	else
		opts = solve_1d_opts(); %defaults
	end;
	
	%Determine whether or not to save the plot
	if nargin > 7
		do_save = 1;
		fn = varargin{3};
	else
		do_save = 0;
	end;
	
	[Sexp, k_quad] = generate_test_data(chi, f, A, ki, zf, noise_ratio, opts);
	
	%Attempt to solve
	[chi_exp,z_exp] = solve_1d(f, Sexp, A, ki, zf, opts);
	
	%Plot solution
	figure;
	hold on;
	plot(z_exp,chi_exp);
	plot(z_exp, arrayfun(chi, z_exp),'r');
	legend('Experimental solution','Exact solution');
	hold off;

	%Optionally, save the plot
	if do_save
		print(fn);
	end;
end
