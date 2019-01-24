clear all
close all

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
maxL=200;
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
    plot([recoveries(i) recoveries(i)],[-0.5 0.5],'--k','LineWidth',1.5)
    ylim([-0.1 0.1])
    if i==4
        xlabel('Time')
    end
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












