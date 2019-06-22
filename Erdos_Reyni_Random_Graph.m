function [ A, G ] = Erdos_Reyni_Random_Graph( n, p )
% Generate Erdos Reyni Random Graphs
%
% Input: 
%       n    numebr of vertices
%       p    probablity of an edge between two vertices
%
% Output:
%       G    an Erdos Reyni Random Graph
%
% Xiaozhe Hu (xiaozhe.hu@tufts.edu), Department of Mathematics, Tufts University
% 02/11/2016

% generate the adjancency matrix
A = rand(n,n) < p;
A = sparse(A);
A = triu(A, 1);
A = A + A';

% generate random graph
G = graph(A);


end

