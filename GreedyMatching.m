function flag = GreedyMatching(d1)
% function flag = GreedyMatching(d1) generates a greedy matching of a graph
% represented by its incidence matrix d1.

% ex: if flag = [1 3 2 2 1], we have 3 matching pairs: (1,5),(3,4),(2)
% More info on matching: https://en.wikipedia.org/wiki/Matching_(graph_theory)#Maximum_matching

% Author: Zian Jiang (zj444@nyu.edu)

[N,~] = size(d1);
flag = zeros(1,N);
k = 0;

for i = 1:N
    if flag(i) == 0
        k = k+1;
        % find nodes that are connected to node i
        neighbor = find(d1(i,:)~=0);
        for j = 1:length(neighbor)
            % should be 1 and -1,(v1,)(v2,)
            curr = find(d1(:,neighbor(j))~=0);
            if flag(curr)==0
                flag(curr)=k;
                break
            end
        end
        flag(i)=k;
    end
end