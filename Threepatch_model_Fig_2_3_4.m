clear all
close all


%% create data
use Threepatch_periodic_add
c=3.1
n=[0.02;0.02;0.02]
eulerstep=0.01;
solver euler eulerstep
out N
series_length=20000;
stabil
data=time(series_length*eulerstep+10);
data=data(:,2:4)-mean(data(:,2:4));
Xs=-1:0.01:1;
Ys=Xs;
Zs=zeros(length(Xs));
ACs=zeros(length(Xs));
Vars=zeros(length(Xs));

for i = 1:length(Xs)
        for j = 1:length(Ys)
            if sqrt(Xs(i)^2+Ys(j)^2)<=1 %only inside circle
                
                Zs(j,i)=sqrt(1-(Xs(i)^2+Ys(j)^2));
                
                PC=[Xs(i) Ys(j) Zs(j,i)];
                y=PC*data';
                                
                ACs(j,i)=corr(y(2:end)',y(1:end-1)');   
                Vars(j,i)=var(y);                            
            end
     end
end

min_AC=min(ACs(ACs>0));
max_AC=max(max(ACs));
ACs(ACs==0)=NaN;

[Wmaf, expl_AC]=MAF(data);
for i=1:length(Wmaf(1,:))
   if Wmaf(3,i)<0
       Wmaf(:,i)=Wmaf(:,i)*-1;
   end   
end

proj_MAF1=data(1:50:end,:)*Wmaf(:,1);
proj_MAF1=proj_MAF1-mean(proj_MAF1);
proj_MAF1=proj_MAF1/std(proj_MAF1);
proj_MAF3=data(1:50:end,:)*Wmaf(:,3);
proj_MAF3=proj_MAF3-mean(proj_MAF3);
proj_MAF3=proj_MAF3/std(proj_MAF3);

autocorr_MAF1=corr(proj_MAF1(2:end),proj_MAF1(1:end-1));
autocorr_MAF3=corr(proj_MAF3(2:end),proj_MAF3(1:end-1));
annotation_AC1=sprintf('  AR(1)= %0.3f',(autocorr_MAF1))
annotation_AC3=sprintf('  AR(1)= %0.3f',(autocorr_MAF3))


%% perform perturbation experiments
use Threepatch_periodic_add
out N
c=3.1
n=0.0

N=[6;6;6];
stabil
Neq=N;
pert_size=0.06;
N_MAF1=N+pert_size*Wmaf(:,1); 
N_MAF2=N+pert_size*Wmaf(:,2);
N_MAF3=N+pert_size*Wmaf(:,3);

S1=time(200);

N=N_MAF1;
S2=time(500);

dist_to_eq=S2(:,2:4)-ones(length(S2(:,1)),1)*Neq';
eucl_dist=sqrt(sum(dist_to_eq.^2,2));
recovery1(1)=min(find(eucl_dist<pert_size*0.1))
recovery1(2)=min(find(eucl_dist<pert_size*0.5))
recovery1(3)=min(find(eucl_dist<pert_size*0.9))

S_MAF1=[S1(end-40:end,2:4); S2(1:360,2:4)];


N=N_MAF2;
S3=time(500);

dist_to_eq=S3(:,2:4)-ones(length(S3(:,1)),1)*Neq';
eucl_dist=sqrt(sum(dist_to_eq.^2,2));
recovery2(1)=min(find(eucl_dist<pert_size*0.1))
recovery2(2)=min(find(eucl_dist<pert_size*0.5))
recovery2(3)=min(find(eucl_dist<pert_size*0.9))

S_MAF2=[S1(end-40:end,2:4); S3(1:360,2:4)];

N=N_MAF3;
S4=time(500);
dist_to_eq=S4(:,2:4)-ones(length(S4(:,1)),1)*Neq';
eucl_dist=sqrt(sum(dist_to_eq.^2,2));

recovery3(1)=min(find(eucl_dist<pert_size*0.1))
recovery3(2)=min(find(eucl_dist<pert_size*0.5))
recovery3(3)=min(find(eucl_dist<pert_size*0.9))

S_MAF3=[S1(end-40:end,2:4); S4(1:360,2:4)];


%% Create figure 2
figure('pos',[100 100 1200 600])
subplot(5,5,[1:5])
plot(data-mean(data),'LineWidth',1.7)
legend('X','Y','Z')
xlim([0 20000])
dim = [0.13 0.636 0.2 0.3];
str = {'A'};
annotation('textbox',dim,'EdgeColor','none','String',str,'FitBoxToText','on','FontSize',14);

a1=subplot(5,5,[6 7 11 12 16 7 21 22])
pos1=get(a1,'Position');
pcolor(Xs,Ys,ACs);
set(gca, 'CLim', [min_AC,max_AC])
pbaspect([1 1 1])
colorbar('westoutside')
hold on;
shading interp; 
title('Autocorrelation')
xlabel('X')
ylabel('Y')
hold on
scatter(0,0,'k','filled')
dim = [0.253 0.155 0.2 0.3];
str = {'Z'};
annotation('textbox',dim,'EdgeColor','none','String',str,'FitBoxToText','on','FontSize',14);
hold on
scatter(Wmaf(1,1),Wmaf(2,1),30,'bx','LineWidth',3)
hold on
scatter(Wmaf(1,3),Wmaf(2,3),'rx','LineWidth',3)
dim = [0.13 0.39 0.2 0.3];
str = {'B'};
annotation('textbox',dim,'EdgeColor','none','String',str,'FitBoxToText','on','FontSize',14);
set(a1,'Position',[0.13 0.10 0.28 0.62])

a2=subplot(5,5,[8 13])
pos2=get(a2,'Position');
scatter(proj_MAF1(1:end-1),proj_MAF1(2:end),15,'k')
hold on
scatter(-1.45,-3.95,30,'bx','LineWidth',3)
hold on
plot([5 5],[-5 5],'k')
hold on
plot([-5 5],[5 5],'k')
hold on
plot([-5 5],[-5 5],'-k')
dim = [0.485 0.18 0.1 0.3];
str = {annotation_AC1};
annotation('textbox',dim,'String',str,'FitBoxToText','on');
pbaspect([1 1 1])
x=xlabel('MAF1(t)')
y=ylabel('MAF1(t+1)')
set(gca,'ytick',[])
set(gca,'xtick',[])
xlim([-5 5])
ylim([-5 5])
dim = [0.445 0.39 0.2 0.3];
str = {'C'};
annotation('textbox',dim,'EdgeColor','none','String',str,'FitBoxToText','on','FontSize',14);
set(a2,'Position',[0.4456 0.42 0.1300 0.28])

a3=subplot(5,5,[18 23])
pos3=get(a3,'Position');
scatter(proj_MAF3(1:end-1),proj_MAF3(2:end),15,'k')
hold on
scatter(-1.45,-3.95,'rx','LineWidth',3)
hold on
plot([5 5],[-5 5],'k')
hold on
plot([-5 5],[5 5],'k')
hold on
plot([-5 5],[-5 5],'-k')
dim = [0.485 0.134 0.2 0.3];
str = {annotation_AC3};
annotation('textbox',dim,'String',str,'FitBoxToText','on','VerticalAlignment', 'bottom');
pbaspect([1 1 1])
xlabel('MAF3(t)')
ylabel('MAF3(t+1)')
set(gca,'ytick',[])
set(gca,'xtick',[])
xlim([-5 5])
ylim([-5 5])
dim = [0.445 0.09 0.2 0.3];
str = {'D'};
annotation('textbox',dim,'EdgeColor','none','String',str,'FitBoxToText','on','FontSize',14);
set(a3,'Position',[0.4456 0.118 0.1300 0.28])

a4=subplot(5,5,[9 10 14 15])
pos4=get(a4,'Position');
plot(S_MAF1-Neq','LineWidth',1.7)
hold on
quiver(40,0,max(recovery1),0,'MaxHeadSize',0.0013,'AutoScale','off','Color','k','LineWidth',2)
hold on
scatter(100,0.36,'bx','LineWidth',3)
ylabel('abundances')
ylim([-.04 0.06])
xlim([0 200])
dim = [0.618 0.39 0.2 0.3];
str = {'E'};
annotation('textbox',dim,'EdgeColor','none','String',str,'FitBoxToText','on','FontSize',14);
set(a4,'Position',[0.6184 0.43 0.2866 0.26])


a5=subplot(5,5,[19 20 24 25])
pos5=get(a5,'Position');
plot(S_MAF3-Neq','LineWidth',1.7)
hold on
quiver(40,0,max(recovery2),0,'MaxHeadSize',0.0025,'AutoScale','off','Color','k','LineWidth',2)
hold on
scatter(100,0.50,'rx','LineWidth',3)
xlabel('time')
ylabel('abundances')
ylim([-.06 0.04])
xlim([0 200])
dim = [0.618 0.09 0.2 0.3];
str = {'F'};
annotation('textbox',dim,'EdgeColor','none','String',str,'FitBoxToText','on','FontSize',14);
set(a5,'Position',[0.6184 0.13 0.2866 0.26])


%% Figure 4A
figure('position', [550, 450, 240, 180])
bar(diag(expl_AC))

%% Figure 3 top row
recoveries=[recovery1;recovery2;recovery3];

figure
subplot(3,1,1)
bar(recoveries(:,1))
subplot(3,1,2)
bar(recoveries(:,2))
subplot(3,1,3)
bar(recoveries(:,3))




