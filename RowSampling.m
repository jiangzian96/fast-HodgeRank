function [As,S,test] = RowSampling(A,sample_size,mode)
% Given A = W*d1', sample rows of A to generate As = S*A
% Sampling modes:
%   1. weighted sampling
%   2. uniform sampling

% Author: Zian Jiang (zj444@nyu.edu)

if mode == "weighted"
    [n_edges,n_nodes] = size(A);
    % PDF
    SqFroNorm = norm(A,'fro')^2;
    total = vecnorm(A,2,2).^2;
    p = total/SqFroNorm;
    % CDF
    cdf = cumsum(p);
    cdf = [0;cdf];
    % sampling
    S = sparse(sample_size,n_edges);
    test = zeros(sample_size,1);
    
     for t = 1:sample_size
        trial = rand;
        result = binary_search(cdf,trial);
        test(t) = result;
     end
     
     for t = 1:sample_size
        S(t,test(t)) = 1/sqrt(sample_size*p(test(t)));
     end
    As = S*A;
end


if mode == "uniform"
    [n_edges,n_nodes] = size(A);
    p = 1/n_edges;
    test = rand(sample_size,1)/p;
    test = ceil(test);
    S = sparse(sample_size,n_edges);
    As = A(test,:);
end