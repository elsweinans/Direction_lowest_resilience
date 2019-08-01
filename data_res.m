function data_res(data,CI,nrsamples,resolutions,blocksize,method)
    nr_vars=size(data,2);
    nr_res=length(resolutions);
    upper=zeros(nr_res,nr_vars);
    lower=zeros(nr_res,nr_vars);
    medians=zeros(nr_res,nr_vars);
    for k=1:nr_res
        res = resolutions(k);
        starts=randi([1 2000001-res*blocksize],1,nrsamples);
        MAFs=zeros(nrsamples,nr_vars);
        for i=1:nrsamples
            data2=data(starts(i):res:starts(i)+(blocksize-1)*res,:);
            if method == 'MAF'
                Wmaf=MAF(data2);
                if Wmaf(1,1)<0
                    Wmaf(:,1)=Wmaf(:,1)*-1;
                end
                MAFs(i,:)=Wmaf(:,1);
            elseif method == 'PCA'
                C=cov(data2);
                [V E]=eig(C);
                d=diag(E);
                idx=find(d==max(d));
                PC=V(:,idx);
                if PC(1)<0
                    PC=PC*-1;
                end
                MAFs(i,:)=PC;
            end
            
        end
        
        for j=1:nr_vars
            A=sort(MAFs(:,j));
            upper(k,j)=A(CI(2));
            lower(k,j)=A(CI(1));
            medians(k,j)=median(A);
        end
        
    end
    
    figure
    hold on
    clrs=['r','g','b','k','c'];
    clrs=clrs(1:nr_vars);
    x=1:length(resolutions);    
    for k=1:nr_vars        
        fill([resolutions flip(resolutions)],[lower(:,k)' flip(upper(:,k)')],clrs(k),'LineStyle','none','facealpha',0.08)
        plot(resolutions,medians(:,k),'Color',clrs(k))
    end   
    ylim([-1 1])
    xlabel('sampling distance')
    if method=='MAF'
        ylabel('MAF')
    elseif method=='PCA'
        ylabel('PC')
    end    
    

end