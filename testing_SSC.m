% Testing script for all deterministic SSC methods on a Watts-Storgatz graph
% Compare r_cgls and r_SSC to see that we get the same ranking from solving
% using either SSC or CGLS (traditional method)

% Author: Zian Jiang (zj444@nyu.edu)

clear
n = 800;
K = round(5*log(n));
[h,L] = WattsStrogatz(n,K,1);
d1 = GenerateD1(L);
[n_nodes,n_edges] = size(d1);
Param.max_it = 3e3;

% change here
Param.solving_mode = 1;
Param.sampling_mode = "weighted";

Param.Tol = 1e-6;
Param.sample_size = round(4*n_nodes*log(n_nodes));
Vals.D1 = d1;
Vals.F = randi([1 5],n_edges,1);
Vals.ww = rand(n_edges,1);
Vals.X0 = zeros(n_nodes,1);
Vals.p = speye(n_nodes,n_nodes);
W = spdiags(Vals.ww(:),0,length(Vals.ww),length(Vals.ww));
A = W.^(1/2)*Vals.D1';
b = W.^(1/2)*Vals.F;
% solve
[x,iter,~,t,~] = SSC(Param,Vals);
X = cgls(A,b,0,1e-6,3e3,false,Vals.X0);

% compare top 10 result
r_cgls = GenerateTop(x,10);
r_SSC = GenerateTop(X,10);