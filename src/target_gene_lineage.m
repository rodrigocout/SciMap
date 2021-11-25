%% very messy/inefficient code
clear

%% data collection (really inefficient...)
targetgenes = {'COL2A1','ACAN','SOX9'};
days = [1,7,14,28,42];
dataset = {};

for n=1:length(targetgenes)
    gene = targetgenes{n};
    datavec = {};
    for k=1:length(days)
        day = days(k);
        load(['D',num2str(day),'_data.mat']);
        T = readtable(['genes_D',num2str(day),'.csv']);
        f_mat = x;
        nmlist = {};
        for j=1:height(T)
            nm = T{j,2};
            nm = nm{1};
            nmlist{end+1} = nm;
        end
        t = find(strcmp(nmlist, gene));
        datavec{k} = f_mat(t,:);
    end
    dataset{n} = datavec;
end

%% 

for k=1:length(days)
    dat1 = dataset{1}{k} > 0;
    dat2 = dataset{2}{k} > 0;
    dat3 = dataset{3}{k} > 0;
    
    v_inp2 = zeros(1,7);
    for j=1:length(dat1)
       a = dat1(j);
       b = dat2(j);
       c = dat3(j);
       if a && b && c
           v_inp2(7) = v_inp2(7)+1;
           v_inp2(5) = v_inp2(5)+1;
           v_inp2(3) = v_inp2(3)+1;
           v_inp2(1) = v_inp2(1)+1;
       elseif a && b
           v_inp2(2) = v_inp2(2)+1;
           v_inp2(3) = v_inp2(3)+1;
           v_inp2(1) = v_inp2(1)+1;
       elseif b && c
           v_inp2(4) = v_inp2(4)+1;
           v_inp2(5) = v_inp2(5)+1;
           v_inp2(3) = v_inp2(3)+1;
       elseif c && a
           v_inp2(6) = v_inp2(6)+1;
           v_inp2(5) = v_inp2(5)+1;
           v_inp2(1) = v_inp2(1)+1;
       elseif a
           v_inp2(1) = v_inp2(1)+1;
       elseif b
           v_inp2(3) = v_inp2(3)+1;
       elseif c
           v_inp2(5) = v_inp2(5)+1;
       end
    end
    A = [3,3,3];
    I = [1, 1, 1, 1];
    A2 = [v_inp2(1), v_inp2(3), v_inp2(5)];
    I2 = [v_inp2(2), v_inp2(6), v_inp2(4), v_inp2(7)];
    figure
    [H,S] = venn(A,I,'ErrMinMode','ChowRodgers','FaceAlpha', 0.6);
    
    %Now label each zone
    for i = 1:3
        text(S.ZoneCentroid(i,1), S.ZoneCentroid(i,2),num2str(targetgenes{i}))
        text(S.ZoneCentroid(i,1), S.ZoneCentroid(i,2)-0.1,['Total: ', num2str(A2(i))])
    end
    for i=4:7
        text(S.ZoneCentroid(i,1), S.ZoneCentroid(i,2), num2str(I2(i-3)))
    end
    set(gca,'visible','off')
    title(['Day ',days(k)]);

%     dat1(find(dat3)) = [];
%     dat2(find(dat3)) = [];
%     dat = zeros(length(dat1),1);
%     for j=1:length(dat1)
%        if dat1(j) > 0 && dat2(j) > 0
%            dat(j) = 3;
%        elseif dat1(j) > 0
%            dat(j) = 1;
%        elseif dat2(j) > 0
%            dat(j) = 2;
%        else
%            dat(j) = 0;
%        end
%     end
%     figure
%     histogram(dat);
%     xticklabels({'None','COL2A1','ACAN','Both'})
%     xticks([0 1 2 3])
%     title('Expression of chondrocyte markers in cells with expression of SOX9');
end


gene = 3;
targetgene = targetgenes{gene};
datavec = dataset{gene};
ns = 0;
figure
subplot(1,5,1);
for k=1:length(datavec)
   subplot(1,5,k);
   histogram(datavec{k}, 'Normalization','probability');
   %hist(datavec{k});
   title(['Expression of ',targetgene,' at day ',num2str(days(k))]);
   xlabel('Expression level (normalized)');
   ylabel('Part of population');
   h = get(gca,'Children');
   disp(sum(datavec{k} > 0))
   
end
set(gcf,'Position',[100 100 1500 200])

% for k=1:length(datavec)
%     ns = ns+length(datavec{k});
% end
% x = zeros(ns,1);
% y = zeros(ns,1);
% ind = 1;
% for k=1:length(datavec)
%     dat = datavec{k};
%     for j=1:length(dat)
%         x(ind) = k;
%         y(ind) = dat(j);
%         ind = ind+1;
%     end
% end
% 
% 
% figure
% plot(x,y,'o');