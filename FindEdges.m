function [G,s,t,w,edges] = FindEdges(d1,r,weights)
d1_t = d1';
n = length(r);
w = zeros(1,n*(n-1)/2);
flag = zeros(1,length(weights));
s = zeros(1,n*(n-1)/2);
t = zeros(1,n*(n-1)/2);
count = 0;
for i = 1:n
    for j = 1:n
        if i ~= j 
            nodes = [r(i),r(j)];
            test = d1_t(:,nodes);
            edge_index = find(sum(abs(test),2)==2);
           
            %flag(edge_index) = 1;
            if flag(edge_index)==0
                count = count + 1;
                s(count) = r(i);
                t(count) = r(j);
                w(count) = weights(edge_index);
                flag(edge_index) = 1;
                edges(count) = edge_index;
            end
            %flag(edge_index) = 1;
        end
    end
end

names = num2str(r);
names = cellstr(names);
%s = unique(s,'stable');
%t = unique(t,'stable');
G = graph(s,t,w);
G = rmnode(G,setdiff(1:max(r),r));
            