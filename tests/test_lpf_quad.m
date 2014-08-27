%Boilerplate
function ret = test_lpf_quad(varargin)
	test_names = find_tests(mfilename());
	tests={};
	for ii = 1:length(test_names) 
		tests{ii} = str2func(test_names{ii}); 
	end;
	ret = run_tests(tests, test_names);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TEST FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [status, msg] = test_sinc_trivial()
	%A low pass filtered sinc should give us a sinc at a lower frequency
	fs = 3;
	Ts = 1/fs;
	t = (-10:Ts:10)';
	tbar = (t+10)/20;
	x = sinc(t);
	x_l = lpf_quad(x, tbar, 10*pi);
	%x_l should go approximately like 0.5*sinc(0.5t).
	%Appears to, except near the ends we have a bit too much attentuation
	%But very good agreement in the centre
	figure;
	hold on;
	plot(t,sinc(0.5*t),'r');
	plot(t,2*x_l);
	hold off;
	fprintf("\n------- test_sinc_trivial--------\nUser Input required:\n\n");
	ans = input('Do the two graphs (more or less) coincide? Type Y for YES or N for NO ','s');
	fprintf("\n");
	if (ans(1) == 'y' || ans(1) == 'Y')
		status = 0;
		msg = 'PASS';
	else
		status = 1;
		msg = 'FAIL';
	end;
end

function [status, msg] = test_sinc_gauss10()
	%This test is to ensure that the filter works correctly for non trivial
	%quadratures
	
	%A low pass filtered sinc should give us a sinc at a lower frequency
	[pts,weights] = generate_quadrature('gauss10',90); %this is equivalent to previous, but we need extra points
	x = sinc(20*pts - 10);
	x_l = lpf_quad(x, pts, 10*pi);
	%x_l should go approximately like 0.5*sinc(0.5t).
	%Appears to, except near the ends we have a bit too much attentuation
	%But very good agreement in the centre
	figure;
	hold on;
	plot(pts,sinc(10*pts-5),'r');
	plot(pts,2*x_l);
	hold off;
	fprintf("\n------- test_sinc_gauss10--------\nUser Input required:\n\n");
	ans = input('Do the two graphs (more or less) coincide? Type Y for YES or N for NO ','s');
	fprintf("\n");
	if (ans(1) == 'y' || ans(1) == 'Y')
		status = 0;
		msg = 'PASS';
	else
		status = 1;
		msg = 'FAIL';
	end;
end

function [status, msg] = test_lowpass_sine()
	%A low frequency sinusoid should pass through the filter without any change
	t = (0:0.01:1)';
	x = sin(2*pi*t);
	x_l = lpf_quad(x, t, 4*pi); %cutoff is double the sinusoid freq
	[status, msg] = assert_eq(x,x_l,0.001); %Graphs are within 0.1%
end

function [status, msg] = test_lowpass_cosine()
	%This test is essentially same as above, but with a cosine instead
	%of a sine
	%This is done to ensure that the filter handles real and imaginary
	%parts correctly in its internal processing (sine has imaginary dft,
	%cosine has real dft)
	
	t = (0:0.01:1)';
	x = cos(2*pi*t);
	x_l = lpf_quad(x, t, 4*pi); %cutoff is double the sinusoid freq
	[status, msg] = assert_eq(x,x_l,0.001); %Graphs are within 0.1%
end
