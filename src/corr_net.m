function [G,A] = corr_net(f_mat,alpha,mode)

if nargin < 4
    mode = 0;
end

% calculate correlation coefficient and corresponding adjacency matrix +
% network
sz = size(f_mat);
c = corrcoef(f_mat');

% mode separates between soft threshold and hard threshold
if mode
    c(isnan(c)) = 0;
    A = abs(c).^alpha;
else
    A = abs(c) > alpha;
end

A = A - diag(diag(A));

G = graph(A);

end