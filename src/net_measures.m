clear
day = 42;
% network measures

load(['Graph_D',num2str(day),'.mat']);

pr = centrality(G,'pagerank');
plot(pr,'o');
[M,I] = maxk(pr,15);
pr_gene_list = nmlist(I);
disp('Genes ranked by pagerank for naive network');
for k=1:length(pr_gene_list)
  nm = pr_gene_list(k);
  disp(nm{1});  
end
t = find(strcmp(nmlist, 'SOX2'));


