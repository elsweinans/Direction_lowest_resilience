clear all
close all

%% Generate data and calculate MAFs
use Chen_model
P=0.35;
n=0.2 
z1=1;
z2=0;
z3=1;
z4=3;
z5=2;
solver euler 0.001

S=time;
data=S(:,2:6);
data=data-mean(data);

[Wmaf expl_AC]=MAF(data)

size_pert=0.6
Zeq=[0;0;0;0;0];
recs105090=zeros(5,3);

%% Perform perturbation experiments
n=0
figure
maxL=300
 
for i = 1:5
    i
    z1=1+size_pert*Wmaf(1,i);
    z2=0+size_pert*Wmaf(2,i);
    z3=1+size_pert*Wmaf(3,i);
    z4=3+size_pert*Wmaf(4,i);
    z5=2+size_pert*Wmaf(5,i);

    S=time(50,'-s');
    data_pert1=S(:,2:6);
    S1=data_pert1-mean(data_pert1);
    
    dist_to_eq=S1-ones(length(S1(:,1)),1)*Zeq';
    eucl_dist=sqrt(sum(dist_to_eq.^2,2));
    recs105090(i,1)=min(find(eucl_dist<size_pert*0.9));
    recs105090(i,2)=min(find(eucl_dist<size_pert*0.5));
    recs105090(i,3)=min(find(eucl_dist<size_pert*0.1));
    
    figure
    subplot(5,1,i)
    plot(S1(1:maxL,:),'LineWidth',1.5)
    hold on
    plot([recs105090(i,1) recs105090(i,1)],[-0.5 0.5],'-k','LineWidth',1.5)
    plot([recs105090(i,2) recs105090(i,2)],[-0.5 0.5],'--k','LineWidth',1.5)
    plot([recs105090(i,3) recs105090(i,3)],[-0.5 0.5],':k','LineWidth',1.5)
    ylim([-0.5 0.5])
end

%% Figure 4C
expl_AC=diag(expl_AC)/sum(diag(expl_AC))
figure('position', [550, 450, 240, 180])
bar(expl_AC)

%% Figure 3 bottom row
figure
for i=1:3
    subplot(3,1,i)
    bar(recs105090(:,4-i))
end

