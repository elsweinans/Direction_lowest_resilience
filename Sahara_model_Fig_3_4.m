clear all
close all
rng('default')

%% Generate data to calculate MAFs
use Sahara_model_4D_noise

B=60; %100

solver euler 0.01
out V
n=0.02*ones(4,1); %0.08*ones(4,1);
S=time(20010);

data=S(1002:end,2:5);

%save('data_sahara')
load('data_sahara')

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
data_suff_length2(data,[6 95],100,'MAF')
data_suff_length2(data,[6 95],100,'PCA')

data_res(data,[6 95],100,[1 2 4 6 10:20:1000],1000,'MAF')
data_res(data,[6 95],100,[1 2 4 6 10:20:1000],1000,'PCA')

%% 1000 perturbation experiments
circ=randn(1000,g_grind.statevars.dim);
dcirc=sqrt(sum(circ.^2,2));
randnrs=circ./dcirc(:,ones(1,g_grind.statevars.dim));
randnrs=randnrs*pert_size;
recs10=zeros(1000,1);
recs50=zeros(1000,1);
recs90=zeros(1000,1);

for i=1:1000
    i
    V=Veq'+randnrs(i,:);
    S=time(200,'-s');
    data_pert1=S(:,2:end);
    S1=data_pert1-Veq';
    
    dist_to_eq=S1;
    eucl_dist=sqrt(sum(dist_to_eq.^2,2));
    recs10(i)=min(find(eucl_dist<pert_size*0.9));
    recs50(i)=min(find(eucl_dist<pert_size*0.5));
    recs90(i)=min(find(eucl_dist<pert_size*0.1));
end

for j=1:2
    j
    if j==1
        V=Veq+pert_size*Wmaf(:,1);
    else
        V=Veq+pert_size*Wmaf(:,end);
    end
    S=time(200,'-s');
    data_pert1=S(:,2:end);
    S1=data_pert1-Veq';
    
    dist_to_eq=S1;
    eucl_dist=sqrt(sum(dist_to_eq.^2,2));
    recs10M(j)=min(find(eucl_dist<pert_size*0.9));
    recs50M(j)=min(find(eucl_dist<pert_size*0.5));
    recs90M(j)=min(find(eucl_dist<pert_size*0.1));
end

figure
subplot(3,1,1)
hold on
histogram(recs10)
plot([recs10M(1) recs10M(1)],[0 600],'-k','LineWidth',1.5)
plot([recs10M(2) recs10M(2)],[0 600],':k','LineWidth',1.5)
ylim([0 600])
subplot(3,1,2)
hold on
histogram(recs50)
plot([recs50M(1) recs50M(1)],[0 200],'-k','LineWidth',1.5)
plot([recs50M(2) recs50M(2)],[0 200],':k','LineWidth',1.5)
ylim([0 200])
ylabel('frequency')
subplot(3,1,3)
hold on
histogram(recs90)
plot([recs90M(1) recs90M(1)],[0 150],'-k','LineWidth',1.5)
plot([recs90M(2) recs90M(2)],[0 150],':k','LineWidth',1.5)
ylim([0 150])
xlabel('recovery times')






