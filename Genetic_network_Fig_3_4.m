clear all
close all

%% Generate data and calculate MAFs
use Chen_model
P=0.35;
n=0.2; 
z1=1;
z2=0;
z3=1;
z4=3;
z5=2;
solver euler 0.001

S=time;
data=S(:,2:6);
data=data-mean(data);

[Wmaf expl_AC]=MAF(data);

size_pert=0.6;
Zeq=[0;0;0;0;0];
recs105090=zeros(5,3);

%% Perform perturbation experiments
n=0;
figure
maxL=300;
figure 
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
    
    subplot(5,1,i)
    plot(S1(1:maxL,:),'LineWidth',1.5)
    hold on
    plot([recs105090(i,1) recs105090(i,1)],[-0.5 0.5],'-k','LineWidth',1.5)
    plot([recs105090(i,2) recs105090(i,2)],[-0.5 0.5],'--k','LineWidth',1.5)
    plot([recs105090(i,3) recs105090(i,3)],[-0.5 0.5],':k','LineWidth',1.5)
    ylim([-0.5 0.5])
end

%% Figure 4C
expl_AC=diag(expl_AC)/sum(diag(expl_AC));
figure('position', [550, 450, 240, 180])
bar(expl_AC)

%% Figure 3 bottom row
figure
for i=1:3
    subplot(3,1,i)
    bar(recs105090(:,4-i))
end


%% Check data usefullness
data=data(5:5:end,:);

resolutions=[1 2 5 10 50 100 500 1000 2000];
figure
 for j=1:length(resolutions)
    j
    res=resolutions(j);
    data2=data(res:res:end,:);
    if resolutions<110        
        L=100:100:length(data2);
    else
        L=10:10:length(data2);
    end
    MAFS=zeros(length(L),5);
    if resolutions < 9
        for i = 1:5:length(L)
            data1=data2(end-L(i)+1:end,:);
            [Wmaf, ~]=MAF(data1);
            if Wmaf(1,1)<0
                Wmaf(:,1)=Wmaf(:,1)*-1;
            end
            MAFS(i,:)=Wmaf(:,1);    
        end
    else
        for i = 1:length(L)
            data1=data2(end-L(i)+1:end,:);
            [Wmaf, ~]=MAF(data1);
            if Wmaf(1,1)<0
                Wmaf(:,1)=Wmaf(:,1)*-1;
            end
            MAFS(i,:)=Wmaf(:,1);    
        end
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
    title(sprintf('time between points = %i',(res)))

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
    PCS=zeros(length(L),5);
    if resolutions < 9
        for i = 1:5:length(L)
            data1=data2(end-L(i)+1:end,:);
            C=cov(data1);
            [V,E]=eig(C);
            PC=V(:,end);
            if PC(1,1)<0
                PC(:,1)=PC(:,1)*-1;
            end
            PCS(i,:)=PC;    
        end
    else
        for i = 1:length(L)
            data1=data2(end-L(i)+1:end,:);
            C=cov(data1);
            [V,E]=eig(C);
            PC=V(:,end);
            if PC(1,1)<0
                PC(:,1)=PC(:,1)*-1;
        end
        PCS(i,:)=PC;    
        end
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
    title(sprintf('time between points = %i',(res)))

 end

