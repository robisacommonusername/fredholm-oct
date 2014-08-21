function ret = test_fastcall_gauss_kernel(varargin)
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

function [status,msg] = test_fastcall_kernel_trivial()
	Q=0;
	alpha = 3;
	z0 = 0;
	npoints = 20;
	f = fastcall_gauss_kernel(Q,alpha,z0,fastcall_opts('n',npoints,'quad_method','trivial'));
	A = ones(npoints,1);
	[Kd,Kdag,pts,k,z,der] = f(A,1,2,3,0);
	k_i = k(pts);
	z_i = z(pts);
	[zz,kk] = meshgrid(z_i,k_i);
	H = gauss_kernel(Q,alpha,z0);
	
	%For trivial quadrature, weights are all 1/npoints
	Kdd = 1/npoints * arrayfun(@(kkk,zzz) H(kkk,zzz)*der(zzz),kk,zz);
	%keyboard();
	[status,msg] = assert_eq(Kd,Kdd);
end