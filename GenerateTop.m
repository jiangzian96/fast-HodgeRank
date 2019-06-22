function I = GenerateTop(r,n)
% function I = GenerateTop(r,n) generates top n candidates (indexes) given
% solution r

% Author: Zian Jiang (zj444@tufts.edu)

[~,I] = sort(r,'descend');
I = I(1:n);
end
