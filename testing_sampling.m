% Testing script for all randomized methods on an Erdos-Reyni graph
% Compare r and r_t to see if the averaged sampled ranking is far from the solution
% on the full graph using CGLS

% Author: Zian Jiang (zj444@nyu.edu)

clear
n_nodes = 5e2;
p = 0.99;
[~,G] = Erdos_Reyni_Random_Graph(n_nodes,p);
L = laplacian(G);
d1 = GenerateD1(L);
[n_nodes,n_edges] = size(d1);
w = exp(randi([-10 10],n_edges,1));
W = spdiags(w,0,n_edges,n_edges);
x0 = zeros(n_nodes,1);
f = randi([1 5],n_edges,1);
x = cgls(W.^(1/2)*d1',W.^(1/2)*f,0,1e-6,3e3,false,x0);
x_t = zeros(n_nodes,20);
% change here
sampling_time = round(25*log(n_nodes));

for i = 1:sampling_time
    % change here
    sampleMode = "uniform";
    sampleSize = round(16*n_nodes*log(n_nodes));
    
    [As,S] = RowSampling(W.^(1/2)*d1',sampleSize,sampleMode); 
    [x_t(:,i),niter(i),flag,t(i)] = conjgrad(As'*As + 1e-8*speye(n_nodes,n_nodes), d1*W*f, x0, 1e-6, 3000);
    fprintf("Finished %d out of %d\n",i,sampling_time)
end

% average sampled results
x_t = mean(x_t,2);

true_sol = GenerateTop(x,n_nodes);
r = GenerateTop(x,10);
r_t = GenerateTop(x_t,10);