%Convert the kernel described by the function H into a frequency domain
%kernel by discretising with n point filon quadrature. n is number of
%quad points required to integrate H(z) by simpsons rule
function [Kd, Kdag, pts, weights] = galerkinise_kernel(H, ki, zf, n, wc)
	%We need to compute \int_0^zf H(k,z) e^iwz dz
	%For every k in ki, and every w in [-wc/2, wc/2]
	Ts = pi/wc;
	N = ceil(1/Ts);
	if mod(N,2) == 0
		N = N+1;
	end;
	zeta = linspace(0,1,N);
	z = unwarp_data('linear',zeta, zf);
	w = 2*pi/N*[0:((N-1)/2), ((1-N)/2):-1]';
	[pts, weights] = filon_quadrature(n, w);
	%weights is z across row, w down column
	[zz,kk] = meshgrid(z,k); %z across row, k down column
	Hvals = arrayfun(H,kk,zz);
	Kd = Hvals*transpose(weights); %weights are complex, avoid '
	Kdag = Kd';
end
