function data_suff_length2(data,CI,nrsamples,method)

    blocks=2000:2000:200000;
    nrblocks=length(blocks);
    nrvars=size(data,2);
    upperlim=zeros(nrblocks,nrvars);
    lowerlim=zeros(nrblocks,nrvars);
    medians=zeros(nrblocks,nrvars);
    for i = 1:nrblocks
        b=blocks(i);
        randnrs=randi([1 length(data)-b],1,nrsamples);
        Wmafs=zeros(nrsamples,nrvars);

        for j=1:nrsamples
            data2=data(randnrs(j):randnrs(j)+b,:);
            if method=='MAF'
                Wmaf=MAF(data2);
                if Wmaf(1,1)<0
                    Wmaf(:,1)=Wmaf(:,1)*-1;
                end
                Wmafs(j,:)=Wmaf(:,1);
            elseif method=='PCA'
                C=cov(data2);
                [V E]=eig(C);
                d=diag(E);
                idx=find(d==max(d));
                PC=V(:,idx);
                if PC(1)<0
                    PC=PC*-1;
                end
                Wmafs(j,:)=PC;
            end                   
        end

        for k=1:nrvars
            A=sort(Wmafs(:,k));
            upperlim(i,k)=A(CI(2));
            lowerlim(i,k)=A(CI(1));
            medians(i,k)=median(A);
        end        
    end

    figure
    hold on
    clrs=['r','g','b','k','c'];
    clrs=clrs(1:nrvars);
    x=1:nrblocks;    
    for k=1:nrvars        
        fill([blocks flip(blocks)],[lowerlim(:,k)' flip(upperlim(:,k)')],clrs(k),'LineStyle','none','facealpha',0.08)
        plot(blocks,medians(:,k),'Color',clrs(k))
    end   
    ylim([-1 1])
    xlabel('time series length')
    if method=='MAF'
        ylabel('MAF')
    elseif method=='PCA'
        ylabel('PC')
    end    
end



