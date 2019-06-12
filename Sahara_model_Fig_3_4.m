clear all
close all
rng('default')

%% Generate data to calculate MAFs
use Sahara_model_4D_noise

B=60; %100

solver euler 0.01
out V
n=0.02*ones(4,1); %0.08*ones(4,1);
S=time(80010);

data=S(1002:end,2:5);
[Wmaf expl_AC]=MAF(data);

if Wmaf(1,1)<0
   Wmaf(:,1)=Wmaf(:,1)*-1;
end
if Wmaf(1,4)<0
   Wmaf(:,4)=Wmaf(:,4)*-1; 
end

%% Perform perturbation experiments
use Sahara_model_4D
out V
B=60;
stabil
Veq=V;
pert_size=0.06; %2
maxL=250;
colors = [1 0.2 0.2; 0 0.5 1; 0  0.6 0.3; 0.5 0.35 0.35];
recs105090=zeros(4,3);

figure
for i = 1:4
    V=Veq-pert_size*Wmaf(:,i);
    simtime 1 100 5000
    S=time(100,'-s');
    data_pert1=S(:,2:5);
    S1=data_pert1-Veq';
    
    dist_to_eq=S1;
    eucl_dist=sqrt(sum(dist_to_eq.^2,2));
    recoveries(i)=min(find(eucl_dist<pert_size/exp(1)));
    recs105090(i,1)=min(find(eucl_dist<pert_size*0.9));
    recs105090(i,2)=min(find(eucl_dist<pert_size*0.5));
    recs105090(i,3)=min(find(eucl_dist<pert_size*0.1));
   
    subplot(4,1,i)
    for j=1:4
        plot(S1(1:maxL,j),'LineWidth',1.5,'Color',colors(j,:))
        hold on
    end
    plot([recs105090(i,1) recs105090(i,1)],[-0.5 0.5],'k','LineWidth',1.5)
    plot([recs105090(i,2) recs105090(i,2)],[-0.5 0.5],'--k','LineWidth',1.5)
    plot([recs105090(i,3) recs105090(i,3)],[-0.5 0.5],':k','LineWidth',1.5)
    
    ylim([-0.1 0.1])
    if i==4
        xlabel('Time')
    end
    ylabel(sprintf('MAF %i',(i)))
end


%% Figure 3 middle row
figure
for i = 1:3
    subplot(3,1,i)
    bar(recs105090(:,4-i))
end

%% Figure 4B
expl_AC=diag(expl_AC)/sum(diag(expl_AC))
figure
bar(expl_AC)
xlabel('MAF numer')
ylabel('explained autocorrelation')

%% Data usefullness
resolutions=[1 2 5 10 50 100 500 1000 2000];
data=data(10:10:end,:);
data=data(end-200000:end,:);
figure
 for j=1:length(resolutions)
     j
    res=resolutions(j);
    data2=data(res:res:end,:);
    %L=100:500:length(data2);
    if resolutions<110        
        L=100:100:length(data2);
    else
        L=10:10:length(data2);
    end
    MAFS=zeros(length(L),4);
    for i = 1:length(L)        
        data1=data2(end-L(i)+1:end,:);
        [Wmaf, ~]=MAF(data1);
        if Wmaf(1,1)<0
            Wmaf(:,1)=Wmaf(:,1)*-1;
        end
        MAFS(i,:)=Wmaf(:,1);    
    end
    
    subplot(3,3,j)
    hold on
    h=plot(L,MAFS,'LineWidth',2);
    xlim([0 max(L)])
    ylim([-1 1])
    xlabel('number of data points')
    ylabel('MAF 1')
    ax = gca;
    ax.XAxis.Exponent = 0;
    title(sprintf('Spacing between points: %i',(res)))
 end

resolutions=[1 2 5 10 50 100 500 1000 2000];
figure
 for j=1:length(resolutions)
    res=resolutions(j);
    data2=data(res:res:end,:);
    if resolutions<110        
        L=100:100:length(data2);
    else
        L=10:10:length(data2);
    end
    PCS=zeros(length(L),4);
    for i = 1:length(L)
        data1=data2(end-L(i)+1:end,:);
        C=cov(data1);
        [V,E]=eig(C);
        PC=V(:,end);
        if PC(1)<0
            PC=PC*-1;
        end
        PCS(i,:)=PC ;  
    end
    
    subplot(3,3,j)
    hold on
    h=plot(L,PCS,'LineWidth',2);
    xlim([min(L) max(L)])
    ylim([-1 1])
    xlabel('number of data points')
    ylabel('PC 1')
    ax = gca;
    ax.XAxis.Exponent = 0;
    title(sprintf('Spacing between points: %i',(res)))

 end









