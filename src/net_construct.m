function [G] = net_construct(f_mat,alpha,day)
% load data
%day = 42;
%alpha = 0.35;
%load(['D',num2str(day),'_data.mat']);
x=f_mat;
% remove genes with low expression
sz = size(f_mat);
rm = [];
ch = [];

for k=1:sz(1)
   if sum(f_mat(k,:)) < 10
       rm(end+1) = k;
   else
       ch(end+1) = sum(f_mat(k,:));
   end
end

f_mat(rm,:) = [];

% calculate correlation coefficient and corresponding adjacency matrix +
% network
sz = size(f_mat);

c = corrcoef(f_mat');

A = abs(c) > alpha;

A = A - diag(diag(A));

rm2 = ~any(A,2);
A( ~any(A,2), : ) = [];  %rows
A( :, ~any(A,1) ) = [];  %columns
G = graph(A);

% indices from initial table to be used
sz = size(x);
ind = 1:sz(1);
ind(rm) = [];
ind(rm2) = [];

% load names
T = readtable(['genes_D',num2str(day),'.csv']);
nmlist = {};
for k=1:height(T)
    nm = T{k,2};
    nm = nm{1};
    nmlist{end+1} = nm;
end

nmlist(rm) = [];
nmlist(rm2) = [];

%save(['Graph_D',num2str(day),'.mat'],'G','nmlist')
end