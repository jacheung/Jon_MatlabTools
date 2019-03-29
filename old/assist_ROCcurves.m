function [ROCsim,ROCx]=assist_ROCcurves(DVsim,DVCR,varargin)
%DVsim = decision variables for y axis; what you are measuring similarity
%TO in ROC analysis (higher DV means that trial is more similar to the mean
%of this variable) 

%DVCR = decision variables for x axis

%varargin = the density between your upper and lower bound (default will be
%.0001

if isempty(varargin)
    spacing = .0001;
else
    spacing = cell2num(varargin(1));
end

    lb=floor(min(horzcat(DVsim,DVCR))); %find lower bound for ranges of criterion
    if lb == 0; %since the floor of zero = zero and crit starts at lb + .0001 we miss capturing those in the P(DV>Crit)
        lb = -1; %set to -1 to increase window of capture
    end
    ub=ceil(max(horzcat(DVsim,DVCR)));%find upper bound for range of criterion

    crit=(lb:spacing:ub); %spacing value for range 
    ROCx=zeros(1,length(crit));
    ROCsim=zeros(1,length(crit));

for c=1:length(crit) %finding P(DV>criterion)
    ROCx(c)=sum(DVCR>crit(c))/numel(DVCR); %should this be numel of DVCR or DVCR+DVsim???
    ROCsim(c)=sum(DVsim>crit(c))/numel(DVsim);
end
