function P_Vc = create_Vc(n_nodes,flag,P)
% function P_Vc = create_Vc(n_nodes,flag,P) creates the coarse grid Vc such
% that every column of output P_Vc represents one coarse subspace

% Author: Zian Jiang (zj444@nyu.edu)

P_Vc = sparse(n_nodes,length(unique(flag)));
for k = 1:max(flag)
    curr = (flag == k);
    curr_Pi = sum(P(:,curr),2);
    P_Vc(:,k) = curr_Pi;
end