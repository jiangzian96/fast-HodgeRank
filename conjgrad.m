function [u, niter, flag,t] = conjgrad(A, b, x0, tol, maxit)
% SOLVECG   Conjugate Gradients method.
%
%    Input parameters: 
%           A : Symmetric, positive definite NxN matrix 
%           f : Right-hand side Nx1 column vector 
%           s : Nx1 start vector (the initial guess)
%         tol : relative residual error tolerance for break
%               condition 
%     maxiter : Maximum number of iterations to perform
%
%    Output parameters:
%           u : Nx1 solution vector
%       niter : Number of iterations performed
%        flag : 1 if convergence criteria specified by TOL could
%               not be fulfilled within the specified maximum
%               number of iterations, 0 otherwise (= iteration
%               successful).

% Author : Andreas Klimke, Universität Stuttgart
% Version: 1.0
% Date   : May 13, 2003
	
u = x0;         % Set u_0 to the start vector s
r = b - A*x0;   % Compute first residuum
p = r;         
rho = r'*r;
niter = 0;     % Init counter for number of iterations
flag = 0;      % Init break flag

% Compute norm of right-hand side to take relative residuum as
% break condition.
normf = norm(b);
if normf < eps  % if the norm is very close to zero, take the
                % absolute residuum instead as break condition
                % ( norm(r) > tol ), since the relative
                % residuum will not work (division by zero).
  warning(['norm(f) is very close to zero, taking absolute residuum' ... 
					 ' as break condition.']);
	normf = 1;
end
tic
for niter = 1:maxit   % Test break condition
	a = A*p;
	alpha = rho/(a'*p);
	u = u + alpha*p;
	r = r - alpha*a;
	rho_new = r'*r;
	p = r + rho_new/rho * p;
	rho = rho_new;
	if (norm(r)/normf < tol)
        t = toc;
        break
    end
	if (niter == maxit)         % if max. number of iterations
		flag = 1;                   % is reached, break.
        t = toc;
		break
	end
end