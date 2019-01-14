% implemented algorithm from https://arxiv.org/pdf/1502.01073.pdf

% Function that decomposes data into a new space where
% the first column is the vector on which the data has the highest
% autocorrelation.

function [Wmaf K]=MAF(data)
    data=data-mean(data);
    Sz=cov(data);
    [U,D]=eig(Sz);
    X=data*U*D^-0.5*U'; 
    dX=X(1:end-1,:)-X(2:end,:);
    Sd=cov(dX);
    [V,K]=eig(Sd);
    Wmaf=U*D^-0.5*U'*V;
    Wmaf=normc(Wmaf);
end