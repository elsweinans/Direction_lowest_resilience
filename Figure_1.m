clear all
close all
rng('default')

%% figure 1A

x=[-9.95:0.1:10];
y=x;
n=length(x);
z=zeros(n,n);
line1=zeros(length(x),1);
line2=zeros(length(x),1);

for i=1:n
    for j=1:n
       
            z(i,j)=x(i)^2+y(j)^2-x(i)*y(j);
            if x(i)==y(j)
                line1(i)=z(i,j);
            elseif x(i)+y(j)<0.01
                i
                line2(i)=z(i,j);
            end

    end
end 

figure
surf(x,y,z')
shading interp
set(gca,'TickLength',[0 0])
set(gca,'XTick',[])
set(gca,'YTick',[])
set(gca,'ZTick',[])
colormap(pink)
hold on
plot3(x,y,line1,'--r','LineWidth',2)
plot3(x,-y,line2,'r','LineWidth',2)

%% Create timeseries for x and y

n=1000;
x=zeros(1,n);
y=zeros(1,n);
ARx=0.9;
ARy=0.1;
        
for i = 2:n
   x(i)=ARx*x(i-1)+normrnd(0,0.2);
   y(i)=ARy*y(i-1)+normrnd(0,0.2) ;   
end

x=x/std(x);
y=y/std(y);

data=[x' y'];
data=data*[cos(pi/4) -sin(pi/4);sin(pi/4) cos(pi/4)];


%% Calculate autocorrelation for every angle
angles=linspace(0,180,201);
ACs=zeros(1,length(angles));

for i=1:length(angles)      
      angle=angles(i);
      PC=[1;tan(pi*angle/180)];
      PC=PC/norm(PC);
      y=PC'*data';
      ACs(i)=corr(y(2:end)',y(1:end-1)'); 
end

max_ACs=find(ACs==max(ACs));
min_ACs=find(ACs==min(ACs));

x=cos(angles/180*pi);
yp=sqrt(1-x.^2);
ym=-sqrt(1-x.^2);

xx=[x;x];
yyp=[yp;yp];
yym=[ym;ym];
zz=zeros(size(xx));

%% Calculate MAFs
[Wmaf K]=MAF(data);

%% Create figure with AC in all directions
figure
hs=surf(xx,yyp,zz,[fliplr(ACs(1,:));fliplr(ACs(1,:))],'LineWidth' ,5,'EdgeColor','interp')
hold on
surf(xx,yym,zz,[ACs(1,:);ACs(1,:)],'LineWidth' ,5,'EdgeColor','interp')
view(2)
pbaspect([1 1 1])
set(gca, 'xticklabel', []);
set(gca, 'yticklabel', []);
xlabel('X')
ylabel('Y')
colorbar
scatter(data(8:5:end,1),data(8:5:end,2),[],'k','filled')
xlim([-1 1])
ylim([-1 1])
plot([-1 1],[0 0],'k','LineWidth',1.2)
plot([0 0],[-1 1],'k','LineWidth',1.1)
plot([-Wmaf(1,1)*2 Wmaf(1,1)*2],[-2*Wmaf(2,1) 2*Wmaf(2,1)],'r','LineWidth',1.5)
plot([-Wmaf(1,2)*2 Wmaf(1,2)*2],[-2*Wmaf(2,2) 2*Wmaf(2,2)],'--r','LineWidth',1.5)
xlabel('X')
ylabel('Y')
title('Autocorrelation in all directions')

%% Calculate and plot projections on MAF1 and MAF2
MAF1=Wmaf(:,1)'*data';
MAF2=Wmaf(:,2)'*data';

figure
subplot(2,1,1)
plot(MAF1)
set(gca, 'yticklabel', []);
set(gca, 'xticklabel', []);
%xlim([1 300])
ylabel('MAF1')
subplot(2,1,2)
plot(MAF2)
set(gca, 'yticklabel', []);
set(gca, 'xticklabel', []);
%xlim([1 300])
%ylim([11.36 11.95])
ylabel('MAF2')


        