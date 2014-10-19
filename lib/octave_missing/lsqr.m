%LSQR implementation for octave, not currently included as standard
% 
%x = lsqr(A,b)
%lsqr(A,b,tol)
%lsqr(A,b,tol,maxit)
%lsqr(A,b,tol,maxit,M)
%lsqr(A,b,tol,maxit,M1,M2)
%lsqr(A,b,tol,maxit,M1,M2,x0)
%[x,flag] = lsqr(A,b,tol,maxit,M1,M2,x0)
%[x,flag,relres] = lsqr(A,b,tol,maxit,M1,M2,x0)
%[x,flag,relres,iter] = lsqr(A,b,tol,maxit,M1,M2,x0)
%[x,flag,relres,iter,resvec] = lsqr(A,b,tol,maxit,M1,M2,x0)
%[x,flag,relres,iter,resvec,lsvec] = lsqr(A,b,tol,maxit,M1,M2,x0)
% 
% Reference: C. C. Paige & M. A. Saunders, "LSQR: an algorithm for
% sparse linear equations and sparse least squares", ACM Trans.
% Math. Software 8 (1982), 43-71.

function [x,flag,relres,iters,residual] = lsqr(A,b, varargin)
	if is_function_handle(A)
		m = length(b);
		n = length(A(b,'transp'));
		fn_handle = 1;
	else
		fn_handle = 0;
		[m,n] = size(A);
	end;
	if nargin > 2
		tol = varargin{1};
	else
		tol = 1e-6; %Compatible with matlab
	end;
	if nargin > 3
		max_iters = varargin{2};
	else
		max_iters = n+1;
	end;
	preconditioned = 0;
	if nargin > 4
		M1 = varargin{3};
		if M1 ~= []
			preconditioned = 1;
			M = M1;
		end;
	end;
	if nargin > 6
		M2 = varargin{4};
	end;
	
	if nargin > 6
		v = varargin{5};
	else
		v = zeros(n,1);
	end;
	%TODO: include preconditioning
	
	%Initialize iteration.
	x = v; 
	beta = norm(b);
	if (abs(beta-eps) < 2*eps) 
		error('Right-hand side must be nonzero')
	end;
	u = b/beta;
	if fn_handle
		r = A(u, 'transp');
	else
		r = (u'*A)';% A'*u, faster to transpose vector twice than matrix;
	end;
	%Todo - work out what all these vars are for, and give them sensible names
	%Initializations are taken from Hansens code.
	alpha = norm(r);
	v = r/alpha; 
	phi_bar = beta;
	rho_bar = alpha;
	w = v;
	c2 = -1; 
	s2 = 0;
	xnorm = 0;
	z = 0;
	
	%These are our vars, and we know what they mean
	iters = 1;
	residual = 0;
	relres = tol+1; %dummy initialization condition
	flag = 0;
	% Perform Lanczos bidiagonalization
	while (iters < max_iters && relres > tol && flag == 0)
		alpha_old = alpha; beta_old = beta;
	
		% Compute A*v - alpha*u.
		if fn_handle
			p = A(v, 'notransp') - alpha*u;
		else
			p = A*v - alpha*u;
		end;
		beta = norm(p); u = p/beta;
	
		% Compute A'*u - beta*v.
		if fn_handle
			r = A(u, 'transp') - beta*v;
		else
			r = (u'*A)' - beta*v; %Double vector transpose faster than matrix transpose
		end;

		alpha = norm(r); v = r/alpha;
	
		% Construct and apply orthogonal transformation.
		rrho = norm([rho_bar,beta]); c1 = rho_bar/rrho;
		s1 = beta/rrho; theta = s1*alpha; rho_bar = -c1*alpha;
		phi = c1*phi_bar; phi_bar = s1*phi_bar;
	
		% Compute solution norm and residual
		%TODO - this appears to be returning incorrect value for overdetermined
		%systems, although the convergence is correct.
		delta = s2*rrho; gamma_bar = -c2*rrho; rhs = phi - delta*z;
		z_bar = rhs/gamma_bar;
		gamma = norm([gamma_bar,theta]);
		c2 = gamma_bar/gamma; s2 = theta/gamma;
		z = rhs/gamma; xnorm = norm([xnorm,z]);
		old_residual = residual;
		residual = abs(phi_bar);
		relres = residual/xnorm;
		
		%Check for stagnation
		if (abs(residual-old_residual) < 2*eps)
			flag = 3;
		end;
		
		% Update the solution.
		x = x + (phi/rrho)*w; w = v - (theta/rrho)*w;
		
		iters = iters + 1;
	end;
	
	%Check if max iterations exceeded
	if iters > max_iters
		flag = 1;
	end;
	
	%Display message if only one output parameter. Done for Matlab compatibility
	if nargout == 1
		if flag == 0
			printf('lsqr converged at iteration %d to a solution with relative residual %f\n',iters,relres);
		elseif flag == 3
			printf('lsqr stagnated after %d iterations\n',iters);
		else
			printf('lsqr did not converge\n');
		end;
	end;
end

%TODO - write tests. Only need octave compatible tests
