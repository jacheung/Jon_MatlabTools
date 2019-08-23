for rec = 1:length(U)
    array = U{rec};
    [objmask]= assist_touchmasks(array);
    mask = objmask.touch;
    kappas = squeeze(array.S_ctk(17,:,:));
    phases = squeeze(array.S_ctk(5,:,:));
    phasemask = nan(size(phases));
    phasemask(phases<0) = 1;
    thetas = squeeze(array.S_ctk(1,:,:));
    thetas(thetas<1) = nan;
    
    
    kappas = kappas.*phasemask.*mask;
    thetas = thetas.*phasemask.*mask;
    
    combined = [kappas(:) thetas(:)];
    keeps = find(~isnan(combined(:,1)));
    combined=combined(keeps,:);
    
    
    
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
    type = {hits,FA,CR};
    for b=1:length(type)
        filttouchIdx = [find(array.S_ctk(9,:,type{b})==1);find(array.S_ctk(12,:,type{b})==1)];
        filtthetas = array.S_ctk(1,:,type{b});
        tmp = [ceil(filttouchIdx/array.t) filtthetas(filttouchIdx)];
        [sorted4, sortedby4 ,~]=binslin(tmp(:,1),tmp(:,2),'equalE',100,.5,100.5);%plotting p(pole go)
        thetatypes{b} = cellfun(@mean,sorted4);
    end
    
    goranges = [floor(min(cell2mat(thetatypes(:,1)))); ceil(max(cell2mat(thetatypes(:,1))))];
    nogoranges =  [min(cell2mat(thetatypes(:,[2 3]))); max(cell2mat(thetatypes(:,[2 3])))];
    nogoranges = [floor(min(nogoranges(1,:))) ceil(max(nogoranges(2,:)))]';

    %sorting kappa:theta pairs and filtering out anything outside of
    %mean+/-2std of each 2 degree bin
    [sorted4, sortedby4 ,bounds]=binslin(combined(:,2),combined(:,1),'equalE',51,1,101);
    ub = cellfun(@mean,sorted4) + 2.*cellfun(@std,sorted4);
    lb = cellfun(@mean,sorted4) - 2.*cellfun(@std,sorted4);    
    filtsorted = cell(length(sorted4),1);
    for k = 1:length(sorted4)
        withinbounds = intersect(find(sorted4{k}>lb(k)),find(sorted4{k}<ub(k)));
        filtsorted{k} = sorted4{k}(withinbounds) ;
    end
    mids = 2:2:100;
    gobins = intersect(find(mids>goranges(1)),find(mids<goranges(2)));
    nogobins = intersect(find(mids>nogoranges(1)),find(mids<nogoranges(2)));
    
    %filtered kappa:theta pairs plotting in figure 40
    allgoskappas = filtsorted{gobins}; allnogoskappas = filtsorted{nogobins};
    [length(allgoskappas) length(allnogoskappas)]
    
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
end    