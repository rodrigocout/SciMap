day = 42;
alpha = 0.35;
load(['D',num2str(day),'_data.mat']);
f_mat = x;
G = net_construct(f_mat,alpha,day);