day = 1;
load(['D',num2str(day),'_data.mat']);
x = x(1:1000,:);
[G,S1] = corr_net(x,0.8,1);
load("ARACNE_net.mat");
S2 = abs(x);

disp(sum(sum(S1~=0)));
disp(sum(sum(S2~=0)))