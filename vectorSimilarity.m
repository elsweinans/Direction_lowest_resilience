function [AngleInRadians,AngleInDegrees,Similarity,pvalue]=vectorSimilarity(I,N)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% this function determines the similarity between two vectors I and N %%%%
%%%% uses PDF in T. Cai, J. Fan, T. Jiang, %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% The Journal of Machine Learning Research 14, 1837 (2013) %%%%%%%%%%%%%%%
%%%% Jelle Lever, 29/03/2018 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% determine number of dimensions
n=length(I);

if n~=length(N)
    error('vector length not equal')
end

%% determine angle between I and N
CosAngle=dot(I,N)/(norm(I)*norm(N));
AngleInRadians=acos(CosAngle);
AngleInDegrees=acosd(CosAngle);

%% stepwise from 0 to radian
NRsteps=1e4;
NULLtoRADIAN=[AngleInRadians/NRsteps:AngleInRadians/NRsteps:AngleInRadians]';

%% PDF from 0 to radian
for STEPNR=1:NRsteps
    testRADIAN=NULLtoRADIAN(STEPNR,1);
    PDF_NULLtoRADIAN(STEPNR,1)=dot((1/(pi.^0.5)).*(gamma(n./2))/(gamma((n-1)./2)),sin(testRADIAN).^(n-2));
end

%% cumulative distribution function
CDF=cumtrapz(NULLtoRADIAN,PDF_NULLtoRADIAN);

%% one-sided p-value
pvalue=CDF(NRsteps,1);

%% similarity measure
Similarity=1-pvalue;