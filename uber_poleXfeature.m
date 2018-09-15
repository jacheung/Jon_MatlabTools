touchCells = touchCell(U,2,.5);
selectedCells = find(touchCells==1);
[rc] = numSubplots(length(selectedCells));

%%
close all;
% figure(89);clf
% figure(88);clf
keep = zeros(6,6);
window = (1:50);
v=1;
for k = 1:length(selectedCells)
    curr = U{selectedCells(k)};
    touchDur = nan(curr.k,50);
    ppolesPro = nan(curr.k,50);
    ppolesRet = nan(curr.k,50);
    velPoles = nan(curr.k,50);

    for b = 1:curr.k 
        touchIdx = [find(curr.S_ctk(9,:,b)==1) find(curr.S_ctk(12,:,b)==1)];
        touchIdxoff = [find(curr.S_ctk(10,:,b)==1) find(curr.S_ctk(13,:,b)==1)];

        if ~isempty(touchIdx)

            velWindow = touchIdx'+(-5:0);
            tmpVelraw = curr.S_ctk(2,velWindow,b);
            
            tmpVel = reshape(tmpVelraw,size(velWindow,1),size(velWindow,2));
            velPoles(b,1:size(tmpVel,1))=nanmean(tmpVel,2)';
            
            
            durdiff = touchIdxoff-touchIdx;
            touchDur(b,1:length(durdiff)) = durdiff; 
            
            tmpThetasraw = curr.S_ctk(v,touchIdx,b);
            tmpPhase = curr.S_ctk(5,touchIdx,b);
            proKeep = tmpPhase<0;
            
            tmpThetas = tmpThetasraw(proKeep);
            tmpThetasret = tmpThetasraw(~proKeep);
            
            ppolesPro(b,1:length(tmpThetas)) = tmpThetas;
            if ~isempty(tmpThetasret)
            ppolesRet(b,1:length(tmpThetasret)) = tmpThetasret;
            end
        end
    end
 %DURATION of 1st touch, THETA, COUNTS, 1st TOUCH PreTouchVelocity, motor position. 
       counts = [sum(~isnan(ppolesPro),2) + sum(~isnan(ppolesRet),2)];
        mat = [touchDur(:,1) nanmean(ppolesPro,2) counts velPoles(:,1) curr.meta.motorPosition'];


        
    % ALL PLOTTING FEATURES
    
     lmmod = fitlm(curr.meta.motorPosition,nanmean(touchDur(:,1),2));
    durcorrVal = abs(sqrt(lmmod.Rsquared.Ordinary));
    figure(87);subplot(rc(1),rc(2),k);
    if durcorrVal>.5
        scatter(curr.meta.motorPosition, nanmean(touchDur(:,1),2),'.b')
    else
        scatter(curr.meta.motorPosition, nanmean(touchDur(:,1),2),'.k')
    end
    set(gca,'xdir','reverse','xlim',[min(curr.meta.motorPosition) max(curr.meta.motorPosition)],'xtick',[],...
        'ylim',[min(touchDur(:,1)) max(touchDur(:,1))])
    
    lmmod = fitlm(curr.meta.motorPosition,nanmean(ppolesPro,2));
    corrVal = abs(sqrt(lmmod.Rsquared.Ordinary));
    figure(88);subplot(rc(1),rc(2),k);
    if corrVal>.5
        scatter(curr.meta.motorPosition, nanmean(ppolesPro,2),'.b')
    else
        scatter(curr.meta.motorPosition, nanmean(ppolesPro,2),'.k')
    end
    set(gca,'xdir','reverse','xlim',[min(curr.meta.motorPosition) max(curr.meta.motorPosition)],'xtick',[],...
        'ylim',[min(nanmean(ppolesPro,2)) max(nanmean(ppolesPro,2))])

    
    counts = [sum(~isnan(ppolesPro),2) + sum(~isnan(ppolesRet),2)];
    lmmodcounts = fitlm(curr.meta.motorPosition,counts);
    countscorrVal = abs(sqrt(lmmodcounts.Rsquared.Ordinary));
    figure(89);subplot(rc(1),rc(2),k);

    if countscorrVal>.5
        scatter(curr.meta.motorPosition, sum(~isnan(ppolesPro),2),'.b')
        hold on; scatter(curr.meta.motorPosition, sum(~isnan(ppolesRet),2),'.c')
    else
        scatter(curr.meta.motorPosition, sum(~isnan(ppolesPro),2),'.k')
        hold on; scatter(curr.meta.motorPosition, sum(~isnan(ppolesRet),2),[],[.5 .5 .5],'.')
    end
    set(gca,'xdir','reverse','xlim',[min(curr.meta.motorPosition) max(curr.meta.motorPosition)],'xtick',[],...
        'ylim',[0 max(sum(~isnan(ppolesPro),2))])

    
    lmmodvel = fitlm(curr.meta.motorPosition,velPoles(:,1));
    velcorrVal = abs(sqrt(lmmodvel.Rsquared.Ordinary));
    figure(90);subplot(rc(1),rc(2),k);

    if velcorrVal>.5
        scatter(curr.meta.motorPosition, velPoles(:,1),'.b')
    else
        scatter(curr.meta.motorPosition, velPoles(:,1),'.k')
    end
    set(gca,'xdir','reverse','xlim',[min(curr.meta.motorPosition) max(curr.meta.motorPosition)],'xtick',[],...
        'ylim',[min(velPoles(:,1)) max(velPoles(:,1))])

    
    corrMat{k} = corr(mat,'rows','complete');
    
    
end


% lab = {'duration','angle','counts','pre touch vel','pole','spks1:50ms'}
% tmp = cell2mat(corrMat');
% tmp2 = tmp(:,5);
% spkcorr = reshape(tmp2,5,length(corrMat));
% 
% figure(81);clf
% for d = 1:size(spkcorr,1)
%     subplot(2,3,d)
%     histogram(spkcorr(d,:),[-1:.1:1],'facecolor','k','facealpha',1,'normalization','probability')
%     set(gca,'ylim',[0 1],'ytick',0:.5:1,'xlim',[-1 1],'xtick',-1:.5:1)
%     title(lab{d})
% end
% suptitle('spike rate correlation')