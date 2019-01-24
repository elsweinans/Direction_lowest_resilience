
%% perturb in all directions to find the direction of longest recovery
use Threepatch_periodic_det

out N
c=3;
stabil
Neq=N;
xs=-1:0.02:1; %xs=-1:0.5:1; % for inaccurate but fast result
ys=-1:0.02:1; %ys=-1:0.5:1; % for inaccurate but fast result
pert_size=0.6;
rectimes=zeros(length(xs),length(ys));

for i=1:length(xs)
    for j=1:length(ys)   
        if norm([xs(i) ys(j)])<1
        
            PC=[xs(i); ys(j);sqrt(1-(xs(i)^2+ys(j)^2))];
            N=Neq+pert_size*PC;
            solver euler 0.01
            simtime 100 1000 ;
            S=time(100,'-s');
            dist_to_eq=S(:,2:end)-Neq';
            eucl_dist=sqrt(sum(dist_to_eq.^2,2));
            rectimes(i,j)=min(find(eucl_dist<pert_size*0.5));    
        end
    end
end

recs=rectimes;
[loc1,loc2]=ind2sub(size(recs),find(recs == max(max(recs))));
longest_rec=[xs(median(loc1)) xs(median(loc2))];
longest_rec=[longest_rec sqrt(1-sum(longest_rec.^2))];

min_recs=min(recs(recs>0));
max_recs=max(max(recs));
recs(recs==0)=NaN;

recs=recs';

%% Initialize noise levels for triangle figure

griddens=57; % griddens=9 % for inaccurate but fast result

stepsize=1/((griddens-1)*2);
X=[];
Y=[];
noises=[];
height=sqrt(1^2-0.5^2);
for i = 1:griddens
    i
    shift=(i-1)*stepsize;
    X=[X shift:2*stepsize:1-shift];
    Y=[Y ones(1,griddens-(i-1))*(i-1)*height/(griddens-1)];
    noises=[noises [1-(2*stepsize)*(i-1):-2*stepsize:0; 0:2*stepsize:1-(2*stepsize)*(i-1); 2*stepsize*(i-1)*ones(1,griddens-(i-1))]];
end

%% Simulate data for all noise levels, calculate MAF and PCA and compare to longest_rec (direction of longest recovery)
use Threepatch_periodic_add
out N
solver euler 0.01
size_noise=0.2;
MAFs=zeros(length(noises),3);
PCAs=zeros(length(noises),3);
MAF1_similarities=zeros(length(noises),1);
PC1_similarities=zeros(length(noises),1);
for i = 1:length(noises)
    i
    noise=noises(:,i);
    n=size_noise*noise;
    N=Neq;
    simtime 1 3000 3000;
    S=time('-s');
    S=S(200:end,2:end);
    S=S-mean(S);
    Wmaf=MAF(S);
    C=cov(S);
    [V,E]=eig(C);
    MAFs(i,:)=Wmaf(:,1);
    PCAs(i,:)=V(:,end); 
    
    [AngleInRadians,AngleInDegrees,Similarity,pvalue]=vectorSimilarity(Wmaf(:,1),longest_rec);
    [AngleInRadians,AngleInDegrees2,Similarity2,pvalue2]=vectorSimilarity(Wmaf(:,1),-longest_rec);
    MAF1_similarities(i)=max(Similarity,Similarity2);

    [AngleInRadians,AngleInDegrees,Similarity,pvalue]=vectorSimilarity(V(:,end),longest_rec);
    [AngleInRadians,AngleInDegrees2,Similarity2,pvalue2]=vectorSimilarity(V(:,end),-longest_rec);
    PC1_similarities(i)=max(Similarity,Similarity2)  ;
    
end

%% Create MAF and PCA triangle plot (Figure 5)

tri = delaunay(X,Y);
triplot(tri,X,Y)
minval=min(min(MAF1_similarities),min(PC1_similarities));
maxval=max(max(MAF1_similarities),max(PC1_similarities));

rico=-0.5*height/0.25;

figure
trisurf(tri,X,Y,MAF1_similarities)
shading interp;
set(gca, 'XLim', [-0.05,1.05])
set(gca, 'YLim', [-0.05,1])
set(gca, 'CLim', [minval,maxval])
hold on

plot3([0 0.02],[0 rico*0.02],[maxval+1 maxval+1],'Color', [0.85 0.85 0.85],'LineWidth',0.8)
plot3([0.125 0.27],[0.25*height 0.25*height+rico*0.145],[maxval+1 maxval+1],'Color', [0.85 0.85 0.85],'LineWidth',0.8)
plot3([0.25 0.52],[0.5*height 0.5*height+rico*0.27],[maxval+1 maxval+1],'Color', [0.85 0.85 0.85],'LineWidth',0.8)
plot3([0.375 0.77],[0.75*height 0.75*height+rico*0.395],[maxval+1 maxval+1],'Color', [0.85 0.85 0.85],'LineWidth',0.8)
plot3([0.5 1.02],[height height+rico*0.52],[maxval+1 maxval+1],'Color', [0.85 0.85 0.85],'LineWidth',0.8)


plot3([0 0.52],[0 -rico*0.52],[maxval+1 maxval+1],'Color', [0.85 0.85 0.85],'LineWidth',0.8)
plot3([0.25 0.645],[0 -rico*0.395],[maxval+1 maxval+1],'Color', [0.85 0.85 0.85],'LineWidth',0.8)
plot3([0.5 0.77],[0 -rico*0.27],[maxval+1 maxval+1],'Color', [0.85 0.85 0.85],'LineWidth',0.8)
plot3([0.75 0.895],[0 -rico*0.145],[maxval+1 maxval+1],'Color', [0.85 0.85 0.85],'LineWidth',0.8)
plot3([1 1.02],[0 -rico*0.02],[maxval+1 maxval+1],'Color', [0.85 0.85 0.85],'LineWidth',0.8)

plot3([-0.03 1],[0 0],[maxval+1 maxval+1],'Color', [0.85 0.85 0.85],'LineWidth',0.8)
plot3([0.085 0.875],[0.25*height 0.25*height],[maxval+1 maxval+1],'Color', [0.85 0.85 0.85],'LineWidth',0.8)
plot3([0.21 0.75],[0.5*height 0.5*height],[maxval+1 maxval+1],'Color', [0.85 0.85 0.85],'LineWidth',0.8)
plot3([0.335 0.625],[0.75*height 0.75*height],[maxval+1 maxval+1],'Color', [0.85 0.85 0.85],'LineWidth',0.8)
plot3([0.46 0.5],[height height],[maxval+1 maxval+1],'Color', [0.85 0.85 0.85],'LineWidth',0.8)

[C,H]=tricontour(tri,X,Y,MAF1_similarities,[0.99 0.99]);
k = numel(H);
for i = 1:k
    A=H(i).Vertices
    plot3(A(:,1),A(:,2),repmat(maxval+1,length(A(:,1)),1),'Color',[0 0 0],'LineWidth',1.5)
end
[C,H]=tricontour(tri,X,Y,MAF1_similarities,[0.95 0.95]);
k = numel(H);
for i = 1:k
    A=H(i).Vertices
    plot3(A(:,1),A(:,2),repmat(maxval+1,length(A(:,1)),1),':k','LineWidth',1.5)
end

grid off
colorbar
pbaspect([1 1 1])
view(2)

figure
trisurf(tri,X,Y,PC1_similarities)
shading interp;
set(gca, 'XLim', [-0.05,1.05])
set(gca, 'YLim', [-0.05,1])
set(gca, 'CLim', [minval,maxval])
hold on
plot3([0 0.02],[0 rico*0.02],[maxval+1 maxval+1],'Color', [0.85 0.85 0.85],'LineWidth',0.8)
plot3([0.125 0.27],[0.25*height 0.25*height+rico*0.145],[maxval+1 maxval+1],'Color', [0.85 0.85 0.85],'LineWidth',0.8)
plot3([0.25 0.52],[0.5*height 0.5*height+rico*0.27],[maxval+1 maxval+1],'Color', [0.85 0.85 0.85],'LineWidth',0.8)
plot3([0.375 0.77],[0.75*height 0.75*height+rico*0.395],[maxval+1 maxval+1],'Color', [0.85 0.85 0.85],'LineWidth',0.8)
plot3([0.5 1.02],[height height+rico*0.52],[maxval+1 maxval+1],'Color', [0.85 0.85 0.85],'LineWidth',0.8)


plot3([0 0.52],[0 -rico*0.52],[maxval+1 maxval+1],'Color', [0.85 0.85 0.85],'LineWidth',0.8)
plot3([0.25 0.645],[0 -rico*0.395],[maxval+1 maxval+1],'Color', [0.85 0.85 0.85],'LineWidth',0.8)
plot3([0.5 0.77],[0 -rico*0.27],[maxval+1 maxval+1],'Color', [0.85 0.85 0.85],'LineWidth',0.8)
plot3([0.75 0.895],[0 -rico*0.145],[maxval+1 maxval+1],'Color', [0.85 0.85 0.85],'LineWidth',0.8)
plot3([1 1.02],[0 -rico*0.02],[maxval+1 maxval+1],'Color', [0.85 0.85 0.85],'LineWidth',0.8)

plot3([-0.03 1],[0 0],[maxval+1 maxval+1],'Color', [0.85 0.85 0.85],'LineWidth',0.8)
plot3([0.085 0.875],[0.25*height 0.25*height],[maxval+1 maxval+1],'Color', [0.85 0.85 0.85],'LineWidth',0.8)
plot3([0.21 0.75],[0.5*height 0.5*height],[maxval+1 maxval+1],'Color', [0.85 0.85 0.85],'LineWidth',0.8)
plot3([0.335 0.625],[0.75*height 0.75*height],[maxval+1 maxval+1],'Color', [0.85 0.85 0.85],'LineWidth',0.8)
plot3([0.46 0.5],[height height],[maxval+1 maxval+1],'Color', [0.85 0.85 0.85],'LineWidth',0.8)

[C,H]=tricontour(tri,X,Y,PC1_similarities,[0.99 0.99]);
k = numel(H);
for i = 1:k
    A=H(i).Vertices
    plot3(A(:,1),A(:,2),repmat(maxval+1,length(A(:,1)),1),'Color',[0 0 0],'LineWidth',1.5)
end
[C,H]=tricontour(tri,X,Y,PC1_similarities,[0.95 0.95]);
k = numel(H);
for i = 1:k
    A=H(i).Vertices
    plot3(A(:,1),A(:,2),repmat(maxval+1,length(A(:,1)),1),':k','LineWidth',1.5)
end

grid off
colorbar
pbaspect([1 1 1])
view(2)


