function [d1] = GenerateD1(L)
% function [d1] = GenerateD1(L) generates boundary operator d1 given graph
% Laplacian L

% Author: Xiaozhe Hu (xiaozhe.hu@tufts.edu)

y = triu(L,1);

% finds the row and column indices and nonzero values of y
[row_y,col_y,~] = find(y);

n_edges = length(row_y);
[n_nodes,~] = size(y);
row_d1 = [row_y; col_y];

col_d1 = [(1:n_edges)'; (1:n_edges)'];

val_d1 = [ones(n_edges,1);-ones(n_edges,1)];

d1 = sparse(row_d1, col_d1, val_d1,n_nodes,n_edges);


end