%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function takes the uber array, masks out periods of touch, and finds
% the kappa associated with the theta angles. Further bins those kappa and
% theta pairs into specific trial types and behavioral outcomes. 
%
% Whisks are PROTRACTION only, thetas < 1 degree are masked out because of
% artifact from sampling grid, and furthermore only kappa values at theta within
% mean+/-2std as there are outliers that contaminate thetakappa pairs. 
%
% Output is 
% - kappaattheta which is a cell array that bins kappatheta pairs by
% trialoutcomes (hit,FA,CR,miss) 
% - tnumskt is cell array of same structure that lists the specific trial
% number associated with each kappatheta pair
% Written 180320 JC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [kappaattheta,tnumskt] = uber_rollmodel(array)

[objmask]= assist_touchmasks(array);
mask = objmask.touch;
kappas = squeeze(array.S_ctk(17,:,:));
phases = squeeze(array.S_ctk(5,:,:));
prophasemask = nan(size(phases));
ret
prophasemask(phases<0) = 1;
thetas = squeeze(array.S_ctk(1,:,:));
thetas(thetas<1) = nan;


kappas = kappas.*prophasemask.*mask;
thetas = thetas.*prophasemask.*mask;
tnums=ceil(intersect(find(~isnan(kappas)),find(~isnan(thetas)))/4000);

%filtering out any trials with nan in kappa/theta
combined = [kappas(:) thetas(:)];
keeps = intersect(find(~isnan(combined(:,1))),find(~isnan(combined(:,2))));
combined=combined(keeps,:);


combined = [combined tnums];


%Variables for classification
nogos =find(array.meta.trialType ==0);
gos = find(array.meta.trialType ==1);
corr = find(array.meta.trialCorrect ==1);
incorr = find(array.meta.trialCorrect == 0);
hits = intersect(gos,corr);
FA = intersect(nogos,incorr);
CR = intersect(nogos,corr);
miss = intersect(gos,incorr);

%finding theta ranges for trialType/trialCorr
type = {hits,FA,CR,miss};
for b=1:length(type)
    filttouchIdx = [find(array.S_ctk(9,:,type{b})==1);find(array.S_ctk(12,:,type{b})==1)];
    filtthetas = array.S_ctk(1,:,type{b});
    tmp = [ceil(filttouchIdx/array.t) filtthetas(filttouchIdx)];
    [sorted, ~ ,~]=binslin(tmp(:,1),tmp(:,2),'equalE',100,.5,100.5);
    thetatypes{b} = cellfun(@nanmean,sorted);
end

hitranges = [floor(min(cell2mat(thetatypes(:,1)))); ceil(max(cell2mat(thetatypes(:,1))))];
FAranges =  [floor(min(cell2mat(thetatypes(:,2)))); ceil(max(cell2mat(thetatypes(:,2))))];
CRranges =  [floor(min(cell2mat(thetatypes(:,3)))); ceil(max(cell2mat(thetatypes(:,3))))];
missranges =  [floor(min(cell2mat(thetatypes(:,4)))); ceil(max(cell2mat(thetatypes(:,4))))];

%sorting kappa:theta pairs and filtering out anything outside of
%mean+/-2std of each 2 degree bin
[sorted4, ~ ,~]=binslin(combined(:,2),combined(:,1),'equalE',51,1,101);
[sortedtmp, ~ ,~]=binslin(combined(:,2),combined(:,[1 3]),'equalE',51,1,101);
ub = cellfun(@mean,sorted4) + 2.*cellfun(@std,sorted4);
lb = cellfun(@mean,sorted4) - 2.*cellfun(@std,sorted4);
filtsorted = cell(length(sorted4),1);
for k = 1:length(sorted4)
    withinbounds = intersect(find(sorted4{k}>lb(k)),find(sorted4{k}<ub(k)));
    filtsorted{k} = sortedtmp{k}(withinbounds,:) ;
end
mids = 2:2:100;


hitbins = intersect(find(mids>hitranges(1)),find(mids<hitranges(2)));
FAbins = intersect(find(mids>FAranges(1)),find(mids<FAranges(2)));
CRbins = intersect(find(mids>CRranges(1)),find(mids<CRranges(2)));
missbins = intersect(find(mids>missranges(1)),find(mids<missranges(2)));


%filtered kappa:theta pairs plotting in figure 40
hitkappas = filtsorted{hitbins};
FAkappas = filtsorted{FAbins};
CRkappas = filtsorted{CRbins};
if ~isempty(missbins) 
misskappas = filtsorted{missbins};
else 
    misskappas = [nan nan];
end

kappaattheta = {hitkappas(:,1),FAkappas(:,1),CRkappas(:,1),misskappas(:,1)};
tnumskt = {hitkappas(:,2),FAkappas(:,2),CRkappas(:,2),misskappas(:,2)};


allgoskappas = [hitkappas ; misskappas ] ;
allnogokappas = [FAkappas; CRkappas];
%% 
    %plotting for raw theta:dkappa
    figure(30);subplot(2,5,rec)
    scatter(combined(:,2), combined(:,1),'.')

    %plotting for mean+/-2std of the ranges (pretty much filter)
    figure(50);subplot(2,5,rec)
    errorbar(1:length(sorted4),cellfun(@mean,sorted4),2.*cellfun(@std,sorted4),'.k')

    %plotting distributions for go/nogo
    figure(40);subplot(2,5,rec)
    histogram(allgoskappas,'normalization','probability','binedges',-.025 :.0025 : .075)
    hold on; histogram(allnogoskappas,'normalization','probability','binedges',-.025 :.0025 : .075)
    set(gca,'xlim',[min([allgoskappas; allnogoskappas])-.025 max([allgoskappas ;allnogoskappas])+.025])
    title(['go:' num2str(goranges(1)) '-' num2str(goranges(2)) ' nogo:' num2str(nogoranges(1)) '-' num2str(nogoranges(2))  ])




