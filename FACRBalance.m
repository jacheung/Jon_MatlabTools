function [FAx,FAy,CRx,CRy1] = FACRBalance(FAx,CRx)
% Function used to fix skewed datasets where we have more of one trial than
% another. Using datasample, we randomly sample from the less counted
% feature in order to repopulate dataset. Datasample only uses whatever
% values are present in the dataset. It doesnt resample from a distribution
% which could be better 

% JC 170804



FAy = ones(size(FAx,1),1)+1;
CRy1 = ones(size(CRx,1),1)+1;

numCR=size(CRx,1);
numFA=size(FAx,1);
numdiff=numFA-numCR;
if numdiff>0
    CRx2=datasample([CRx CRy1],numdiff);
    CRx = [CRx;CRx2(:,1:end-1)];
    CRy1 = [CRy1;CRx2(:,end)];
else
    FAx2 = datasample([FAx FAy],abs(numdiff));
    FAx = [FAx;FAx2(:,1:end-1)];
    FAy = [FAy;FAx2(:,end)];
end