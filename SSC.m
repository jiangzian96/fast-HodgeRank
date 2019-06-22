function [x,iter,res,t,S,err] = SSC(Param,Vals)
% Successive Subspace Correction solver

% Solving modes:
%   1. Canonical decomposition (stable)
%   2. Graph Matching decomposition (stable)
%   3. Canonical decomposition + coarse grid (stable)
%   4. Canonical decomposition preconditioned by sampled L (stable)
%   5. Canonical decomposition preconditioned by sampled L + coarse grid
%   from sampled L (VERY UNSTABLE)

% Sampling modes:
%   1. weighted sampling: sampling edges out considering edge weights w
%   2. uniform sampling: set p = 1/E

% Param:
%   max iteration
%   solving mode
%   sampling mode
%   tol
%   samplesize
% Val:
%   d1: size VxE
%   f: size E
%   w: size E
%   x0: size V
%   P: each column is length E, representing 1 subspace

% Author: Zian Jiang (zj444@nyu.edu)

maxit = Param.max_it;
solveMode = Param.solving_mode;
sampleMode = Param.sampling_mode;
tol = Param.Tol;
sampleSize = Param.sample_size;
d1 = Vals.D1;
f = Vals.F;
w = Vals.ww;
x0 = Vals.X0;
P = Vals.p;
W = spdiags(w(:),0,length(w),length(w));
[n_nodes,n_edges] = size(d1);
A = W.^(1/2)*d1';
b = W.^(1/2)*f;
x = x0;
res = b - A*x;
err = zeros(maxit+1,1);
err(1) = norm(A'*(b - A*x))/norm(A'*f);
if solveMode == 4||solveMode == 5
   [As,S] = RowSampling(A,sampleSize,sampleMode); 
end
if solveMode == 2|| solveMode == 3
    flag = GreedyMatching(A');
    if solveMode == 3
        P_Vc = create_Vc(n_nodes,flag,P);
    end
end
if solveMode == 5
    flag = GreedyMatching(As');
end
tic
S = 0;
for iter = 1:maxit
    if solveMode == 1 || solveMode == 3
        B = tril(A'*A);
        x = x + B\(A'*(b - A*x));
        res = b - A*x;
        err(iter+1) = norm(A'*(b - A*x))/norm(A'*f);
    end
    if solveMode == 2
        for i = 1:max(flag)
            curr = flag == i;
            B1 = A*P(:,curr);
            B2 = (b - A*x);
            alpha = (B1'*B1)\(B1'*B2);
            x = x + P(:,curr)*alpha;
            res = b - A*x;
            err(iter+1) = norm(A'*(b - A*x))/norm(A'*f);
        end
    end
    if solveMode == 3
        nc = size(P_Vc,2);
        B3 = A*P_Vc;
        alpha = (B3'*B3 + 10^(-6)*speye(nc,nc))\(B3'*(b - A*x));
        x = x + P_Vc*alpha;
        res = b - A*x;
        err(iter+1) = norm(A'*(b - A*x))/norm(A'*f);
    end
    if solveMode == 4
        B = tril(As'*As);
        B = B + 1e-6*speye(n_nodes,n_nodes);
        x = x + B\(A'*(b - A*x));
        res = b - A*x;
        err(iter+1) = norm(A'*(b - A*x))/norm(A'*f);
    end
    if solveMode == 5
        B = tril(As'*As);
        B = B + 1e-6*speye(n_nodes,n_nodes);
        x = x + B\(A'*(b - A*x));
        res = b - A*x;
        P_Vc = sparse(n_nodes,length(unique(flag)));
        for k = 1:max(flag)
            curr = (flag == k);
            curr_Pi = sum(P(:,curr),2);
            P_Vc(:,k) = curr_Pi;
        end
        nc = size(P_Vc,2);
        B3 = A*P_Vc;
        alpha = (B3'*B3 + 10^(-6)*speye(nc,nc))\(B3'*(b - A*x));
        x = x + P_Vc*alpha;
        res = b - A*x;
        err(iter+1) = norm(A'*(b - A*x))/norm(A'*f);
    end
    
    if norm(A'*(b - A*x))/norm(A'*f) < tol
        t = toc;
        err = err(1:iter+1);
        break;
    end
    
end




