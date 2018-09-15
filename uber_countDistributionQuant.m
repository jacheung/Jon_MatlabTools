%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script quantifies the touch count distribution by finding the
% distribution with the lowest RMSE value that matches. Template
% distributions are flat,linear, sigmoids, or left/right shifted sigmoids
%
% Main output are figures and the sensory bias value "acqbias"
% acqbias 0 = flat distribution : 1 = left shifted sigmoid towards go

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



type = {SM,BV};

for d = 2
    U=type{d};
    close all
    for rec = 1:length(U)
        
        [prelickTouches] = assist_predecisionVar(U{rec});
        cpoles = cellfun(@numel,prelickTouches)';
        poles = U{rec}.meta.motorPosition;
        poleranges = U{rec}.meta.ranges;
        gos = find(U{rec}.meta.trialType==1);
        nogos = find(U{rec}.meta.trialType==0);
        [sortedCounts]= binslin(poles,cpoles,'equalE',11,poleranges(1),poleranges(end));
        meantouches = cellfun(@mean,sortedCounts);
        normTouches = meantouches./sum(meantouches);
            nogotmp = sum(normTouches(1:6));
            gotmp = sum(normTouches(6:11));
            acqbiascalc(rec) = 1-(nogotmp./gotmp);
        polepositions = linspace(poleranges(1),poleranges(end),length(sortedCounts));
        
        linearvars =  normalize_var(linspace(0,1,length(meantouches)),min(meantouches),max(meantouches));
        sigmoidvars =  normalize_var(sigmoid(-5:5),min(meantouches),max(meantouches));
        lshiftsigmoidvars = [repmat(sigmoidvars(1),1,3) sigmoidvars(1:length(sigmoidvars)-3)];
        rshiftsigmoidvars = [sigmoidvars(4:length(sigmoidvars)) repmat(sigmoidvars(end),1,3)];
        
        linearrmse=mean((meantouches-linearvars').^2);
        linearslope = linearvars(5)-linearvars(4);
        sigmoidrmse=mean((meantouches-sigmoidvars').^2);
        lsigmoidrmse = mean((meantouches-lshiftsigmoidvars').^2);
        rsigmoidrmse = mean((meantouches-rshiftsigmoidvars').^2);
        
        fitsallRMSE(rec,:) = [linearrmse sigmoidrmse lsigmoidrmse rsigmoidrmse];
        fitslope(rec) = linearslope;
        [~,rmsebestidx] = min(fitsallRMSE,[],2);
        
        
        %% ALL Plotting features
        %POLE X TOUCH COUNTS
        figure(801);subplot(2,5,rec)
        bar(polepositions(1:5),meantouches(1:5),'r');hold on; bar(polepositions(6:10),meantouches(6:10),'b')
        scatter(poles(gos),cpoles(gos),'b');hold on; scatter(poles(nogos),cpoles(nogos),'r')
        hold on; plot(linspace(poleranges(1),poleranges(end),length(meantouches)),meantouches,'k','linewidth',5)
        set(gca,'xdir','reverse','xlim',[poleranges(1) poleranges(end)],'xticklabel',[],'ylim',[-1 15],'ytick',0:5:15)
        
        %REAL TOUCH COUNT DISTRIBUTION x BEST FIT
        figure(8500);subplot(2,5,rec)
        plot(1:length(meantouches),meantouches,'k','linewidth',3)
        if rmsebestidx(rec) == 1
            hold on; plot(1:length(meantouches),sigmoidvars,'r-.')
            if fitslope(rec)>=.5
                hold on; plot(1:length(meantouches),linearvars,'g-.')
                hold on; plot(1:length(meantouches),ones(length(meantouches),1).*mean(meantouches),'-.r')
            else
                hold on; plot(1:length(meantouches),linearvars,'r-.')
                hold on; plot(1:length(meantouches),ones(length(meantouches),1).*mean(meantouches),'-.g')
            end
        elseif rmsebestidx(rec)>=2
            if rmsebestidx(rec) == 2
                    hold on; plot(1:length(meantouches),sigmoidvars,'g-.')
                    hold on; plot(1:length(meantouches),linearvars,'r-.')
                    hold on; plot(1:length(meantouches),ones(length(meantouches),1).*mean(meantouches),'-.r')
            elseif rmsebestidx(rec) == 3
                    hold on; plot(1:length(meantouches),lshiftsigmoidvars,'g-.')
                    hold on; plot(1:length(meantouches),linearvars,'r-.')
                    hold on; plot(1:length(meantouches),ones(length(meantouches),1).*mean(meantouches),'-.r')
            elseif rmsebestidx(rec) == 4
                    hold on; plot(1:length(meantouches),rshiftsigmoidvars,'g-.')
                    hold on; plot(1:length(meantouches),linearvars,'r-.')
                    hold on; plot(1:length(meantouches),ones(length(meantouches),1).*mean(meantouches),'-.r')
                    acqbias(rec) = .5;
            end
        end
        axis square
        set(gca,'xdir','reverse','xlim',[1 length(meantouches)],'xticklabel',[],'ylim',[0 15],'ytick',0:5:15)
        
        %EXAMPLE OF ALL FITS        
        figure(38000);clf
        subplot(1,5,5);plot(1:length(meantouches),lshiftsigmoidvars,'k','linewidth',3)
        set(gca,'xdir','reverse','xlim',[1 length(meantouches)],'xticklabel',[],'ylim',[0 10],'yticklabel',[])
        subplot(1,5,4);plot(1:length(meantouches),sigmoidvars,'k','linewidth',3)
        set(gca,'xdir','reverse','xlim',[1 length(meantouches)],'xticklabel',[],'ylim',[0 10],'yticklabel',[])
        subplot(1,5,3); plot(1:length(meantouches),rshiftsigmoidvars,'k','linewidth',3)
        set(gca,'xdir','reverse','xlim',[1 length(meantouches)],'xticklabel',[],'ylim',[0 10],'yticklabel',[])
        subplot(1,5,2);plot(1:length(meantouches),linearvars,'k','linewidth',3)
        set(gca,'xdir','reverse','xlim',[1 length(meantouches)],'xticklabel',[],'ylim',[0 10],'yticklabel',[])
        subplot(1,5,1); plot(1:length(meantouches),repmat(mean(meantouches),1,length(meantouches)),'k','linewidth',3)
        set(gca,'xdir','reverse','xlim',[1 length(meantouches)],'xticklabel',[],'ylim',[0 10],'yticklabel',[])
        
        
        lshifttmp = lshiftsigmoidvars./sum(lshiftsigmoidvars);
        rshifttmp = rshiftsigmoidvars./sum(rshiftsigmoidvars);
        sigmoidtmp = sigmoidvars./sum(sigmoidvars);
        lintmp = linearvars./sum(linearvars);
 
        fulltmp=[lshifttmp;sigmoidtmp;rshifttmp;lintmp];
        nogotmps = sum(fulltmp(:,1:6),2);
        gotmps = sum(fulltmp(:,6:11),2);
        1-(nogotmps./gotmps)
        
    end
end