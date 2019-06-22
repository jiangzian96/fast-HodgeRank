function result = binary_search(cdf,trial)
% Using binary search to find the correct edge that was sampled

% Author: Zian Jiang (zj444@nyu.edu)

N = length(cdf);
n_edges = N-1;

b = numel(cdf);
a = 0;

n = b;
while true
    p = a + floor(n/2);
    if cdf(p) > trial
        b = p;
    else
        a = p;
    end
    
    n = b-a;
    if n == 1
        result = a;
        break
    end
end
    
end
