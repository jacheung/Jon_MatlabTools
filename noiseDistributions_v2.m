%% Counts distribution
type = {SM,BV};

for d = 2
    colors = {'DarkMagenta','DarkTurquoise'};
    U=type{d};
    close all
    for rec = 1:length(U)
        
        [prelickTouches] = assist_predecisionVar(U{rec});
        cpoles = cellfun(@numel,prelickTouches)';
        ppoles = nan(U{rec}.k,2);
        for b = 1:U{rec}.k
            %FOR THETAS
            touchIdx = [find(U{rec}.S_ctk(9,:,b)==1) find(U{rec}.S_ctk(12,:,b)==1)];
            phaseattouch = U{rec}.S_ctk(5,touchIdx,b);
            touchIdx = touchIdx(phaseattouch<=0);
            
            
            if ~isempty(touchIdx)
                tmpThetas = U{rec}.S_ctk(1,touchIdx,b);
                ppoles(b,:) = [median(tmpThetas) std(tmpThetas)];
                
            end
        end
        
        %find trials with lick
        hits = double(U{rec}.meta.trialCorrect == 1) .* double(U{rec}.meta.trialType ==1);
        miss = double(U{rec}.meta.trialCorrect == 0) .* double(U{rec}.meta.trialType ==1);
        FA = double(U{rec}.meta.trialCorrect == 0) .* double(U{rec}.meta.trialType ==0);
        CR = double(U{rec}.meta.trialCorrect == 1) .* double(U{rec}.meta.trialType ==0);
        lix = hits + FA;
        
        
        %polynomial or linear fits to data
        
        %theta
        polyinputs = sortrows([U{rec}.meta.motorPosition'  ppoles(:,1)]);
        polyinputs(isnan(polyinputs(:,2)),:)=[];
        [coefftheta, ~ , mutheta] = polyfit(polyinputs(:,1),polyinputs(:,2),2);
        thetafillscomplete = polyval(coefftheta,U{rec}.meta.ranges(1):500:U{rec}.meta.ranges(2),[],mutheta);
        
        %thetastd
        polyinputsstd = sortrows([U{rec}.meta.motorPosition'  ppoles(:,2)]);
        ogstd = polyinputsstd(:,2);
        polyinputsstd(isnan(polyinputsstd(:,2)),:)=[];
        [coeff, ~ , mu] = polyfit(polyinputsstd(:,1),polyinputsstd(:,2),1);
        stdfits = polyval(coeff,sort(U{rec}.meta.motorPosition),[],mu);
        
        %counts
        polyinputscounts = sortrows([U{rec}.meta.motorPosition' cpoles(:,1)]);
        [coeffcounts, ~ , mucounts] = polyfit(polyinputscounts(:,1),polyinputscounts(:,2),1);
        countssigtest = polyval(coeffcounts,sort(U{rec}.meta.motorPosition),[],mucounts);
        countsfits = polyval(coeffcounts,U{rec}.meta.ranges(1):500:U{rec}.meta.ranges(2),[],mucounts);
        
        %licks
        polyinputslix = sortrows([U{rec}.meta.motorPosition' lix']);
        [coefflix, ~ , mulix] = polyfit(polyinputslix(:,1),polyinputslix(:,2),1);
        lixsigtest = polyval(coefflix,sort(U{rec}.meta.motorPosition),[],mulix);
        lixfits = polyval(coefflix,U{rec}.meta.ranges(1):500:U{rec}.meta.ranges(2),[],mulix);
        
        
        poles = U{rec}.meta.motorPosition;
        [~,idx] = sort(poles);
        poleranges = U{rec}.meta.ranges;
        gos = find(U{rec}.meta.trialType==1);
        nogos = find(U{rec}.meta.trialType==0);
        
        
        % Pole position x counts
        modelvals = fitlm(countssigtest,cpoles(idx));
        pvalcounts = modelvals.Coefficients{2,4};
        countsrsq= modelvals.Rsquared.ordinary;
        
        [sorted]= binslin(poles,cpoles,'equalE',12,poleranges(1),poleranges(end));
        [sortedlix]= binslin(poles,lix','equalE',12,poleranges(1),poleranges(end));
        
        meanlix = cellfun(@mean,sortedlix);
        meantouches = cellfun(@mean,sorted);
        %         figure(803);clf;
        %         scatter(meantouches,meanlix);
        %         set(gca,'xlim',[0 15],'ytick',0:.5:1)
        %
        
        stdtouches = cellfun(@std,sorted);
        polepositions = linspace(poleranges(1),poleranges(end),length(sorted));
        onestdcounts = mean(stdtouches);
        polebins = U{rec}.meta.ranges(1):500:U{rec}.meta.ranges(2);
        
        %POLE X TOUCH COUNTS
        figure(801);subplot(2,5,rec)
        bar(polepositions(1:5),meantouches(1:5),'r');hold on; bar(polepositions(6:10),meantouches(6:10),'b')
        scatter(poles(gos),cpoles(gos),'b');hold on; scatter(poles(nogos),cpoles(nogos),'r')
        hold on; plot(linspace(poleranges(1),poleranges(end),length(meantouches)),meantouches,'k','linewidth',5)
        %         hold on; plot(U{rec}.meta.ranges(1):500:U{rec}.meta.ranges(2),countsfits,'k','linewidth',3)
%
        set(gca,'xdir','reverse','xlim',[U{rec}.meta.ranges(1) U{rec}.meta.ranges(end)],'xticklabel',[],'ylim',[-1 15],'ytick',0:5:15)
        

        
        
        linearvars =  normalize_var(linspace(0,1,length(meantouches)),min(meantouches),max(meantouches));
        linearSSR = sum((linearvars-mean(meantouches)).^2);
        linearSSE = sum((meantouches - linearvars').^2);
        linearSST = linearSSR+linearSSE;
        linearrsq = 1-linearSSE/linearSST;
        linearrmse=mean((meantouches-linearvars').^2);
        linearslope = linearvars(5)-linearvars(4);
        
        sigmoidvars =  normalize_var(sigmoid(-5:5),min(meantouches),max(meantouches));
        sigmoidSSR = sum((sigmoidvars-mean(meantouches)).^2);
        sigmoidSSE = sum((meantouches - sigmoidvars').^2);
        sigmoidSST = sigmoidSSR+sigmoidSSE;
        sigmoidrsq = 1-sigmoidSSE/sigmoidSST;
        sigmoidrmse=mean((meantouches-sigmoidvars').^2);
        
        lshiftsigmoidvars = [repmat(sigmoidvars(1),1,3) sigmoidvars(1:length(sigmoidvars)-3)];
        lsigmoidrmse = mean((meantouches-lshiftsigmoidvars').^2);
        rshiftsigmoidvars = [sigmoidvars(4:length(sigmoidvars)) repmat(sigmoidvars(end),1,3)];
        rsigmoidrmse = mean((meantouches-rshiftsigmoidvars').^2);
        
        linearfit = fitlm(meantouches,normalize_var(linspace(0,1,length(meantouches)),min(meantouches),max(meantouches)));
        sigmoidfit = fitlm(meantouches,normalize_var(sigmoid(-5:5),min(meantouches),max(meantouches)));
        
        fitsallRSQ(rec,:) = [linearfit.Rsquared.ordinary sigmoidfit.Rsquared.ordinary];
        fitsallRMSE(rec,:) = [linearrmse sigmoidrmse lsigmoidrmse rsigmoidrmse];
        fitslope(rec) = linearslope;
        [~,rmsebestidx] = min(fitsallRMSE,[],2);
        
        %DEFINING TOUCH COUNT DISTRIBUTIONS
        figure(8500);subplot(2,5,rec)
        plot(1:length(meantouches),meantouches,'k','linewidth',3) 
        if rmsebestidx(rec) == 1
            hold on; plot(1:length(meantouches),sigmoidvars,'r-.')
            if fitslope(rec)>=.5
                hold on; plot(1:length(meantouches),linearvars,'g-.')
                hold on; plot(1:length(meantouches),ones(length(meantouches),1).*mean(meantouches),'-.r')
                acqbias(rec) = .25;
            else
                hold on; plot(1:length(meantouches),linearvars,'r-.')
                hold on; plot(1:length(meantouches),ones(length(meantouches),1).*mean(meantouches),'-.g')
                acqbias(rec) = 0;
            end
        elseif rmsebestidx(rec)>=2
            if rmsebestidx(rec) == 2
            hold on; plot(1:length(meantouches),sigmoidvars,'g-.')
            hold on; plot(1:length(meantouches),linearvars,'r-.')
            hold on; plot(1:length(meantouches),ones(length(meantouches),1).*mean(meantouches),'-.r')
            acqbias(rec) = .75;
            elseif rmsebestidx(rec) == 3
            hold on; plot(1:length(meantouches),lshiftsigmoidvars,'g-.')
            hold on; plot(1:length(meantouches),linearvars,'r-.')
            hold on; plot(1:length(meantouches),ones(length(meantouches),1).*mean(meantouches),'-.r')
            acqbias(rec) = 1;
            elseif rmsebestidx(rec) == 4 
            hold on; plot(1:length(meantouches),rshiftsigmoidvars,'g-.')
            hold on; plot(1:length(meantouches),linearvars,'r-.')
            hold on; plot(1:length(meantouches),ones(length(meantouches),1).*mean(meantouches),'-.r')
            acqbias(rec) = .5;
            end
        end
  
        set(gca,'xdir','reverse','xlim',[1 length(meantouches)],'xticklabel',[],'ylim',[0 15],'ytick',0:5:15)
        

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
        
        
        
        %         [fitresult, gof] = createFitstest(poles, cpoles);
        %
        %         gofvals(rec,:) = [gof(:).rsquare];
        
        %         pause
        %         %POLE X TOUCH COUNTS STD
        %         [coeffcountsstd, ~ , mucountsstd] = polyfit(polepositions',stdtouches,1);
        %         stdcountsfits = polyval(coeffcountsstd,polepositions,[],mucountsstd);
        %         modelvals = fitlm( stdcountsfits ,stdtouches);
        %         pvalcountsstd = modelvals.Coefficients{2,4};
        %
        %         subplot(2,2,3)
        %         scatter(polepositions(1:25),stdtouches(1:25),'r')
        %         hold on; scatter(polepositions(26:50),stdtouches(26:50),'b')
        %         hold on; plot(polepositions,stdcountsfits,'k');
        %         set(gca,'xdir','reverse','xlim',[U{rec}.meta.ranges(1) U{rec}.meta.ranges(end)],'xticklabel',[],'ylim',[-1 10],'ytick',0:5:15)
        %         ylabel('std touches')
        %         title(['pval=' num2str(pvalcountsstd)  '          mean std = ' num2str(round(nanmean(stdtouches),2))])
        %
        %
        %         %POLE X THETA
        %         subplot(2,2,2)
        %         errorbar(poles(gos),ppoles(gos,1),ppoles(gos,2),'o')
        %         hold on; errorbar(poles(nogos),ppoles(nogos,1),ppoles(nogos,2),'ro')
        %         hold on; plot(U{rec}.meta.ranges(1):500:U{rec}.meta.ranges(2),thetafillscomplete,'k','linewidth',3)
        %         ylabel('theta')
        %         set(gca,'xdir','reverse','xticklabel',[],'xlim',[U{rec}.meta.ranges(1) U{rec}.meta.ranges(end)])
        %
        %         %checking to see how much STD of theta varries across pole positions
        %         modelvals = fitlm(stdfits,ogstd);
        %         pval = modelvals.Coefficients{2,4};
        %         ordrsq = modelvals.Rsquared.Ordinary;
        %         onestdtheta = round(mean(stdfits),2);
        %
        %         subplot(2,2,4)
        %         scatter(poles(gos),ppoles(gos,2),'b')
        %         hold on; scatter(poles(nogos),ppoles(nogos,2),'r')
        %         hold on; plot(sort(poles),stdfits,'k')
        %         set(gca,'xdir','reverse','xticklabel',[],'xlim',[U{rec}.meta.ranges(1) U{rec}.meta.ranges(end)])
        %         ylabel('std theta')
        %         title(['pval=' num2str(pval) '     mean std=' num2str(onestdtheta)])
        %
        
    end
    
    %% generating pseudo values
    
    
    
    
    motorbounds = [U{rec}.meta.ranges(1) mean(U{rec}.meta.ranges) U{rec}.meta.ranges(2)];
    optboundstheta = polyval(coefftheta,motorbounds,[],mutheta);
    optboundscounts = -.5:1:15.5;
    
    sigvaluestheta = [0 onestdtheta/2 onestdtheta onestdtheta*1.5 onestdtheta*2 onestdtheta*2.5 onestdtheta*3];
    sigvaluescounts = [0 onestdcounts/2 onestdcounts onestdcounts*1.5 onestdcounts*2 onestdcounts*2.5 onestdcounts*3];
    
    figure(380);clf
    figure(381);clf
    for p=1:length(sigvaluestheta)
        sampmotornogo = datasample(motorbounds(1):250:motorbounds(2),5000);
        sampthetanogo = polyval(coefftheta,sampmotornogo,[],mutheta);
        sampcountsnogo = polyval(coeffcounts,sampmotornogo,[],mucounts);
        for g = 1:length(sampmotornogo)
            nogodistthetadraw = makedist('normal','mu',sampthetanogo(g),'sigma',sigvaluestheta(p));
            nogodistcountsdraw = makedist('normal','mu',sampcountsnogo(g),'sigma',sigvaluescounts(p));
            nogodisttheta(g) = random(nogodistthetadraw,1,1);
            nogodistcounts(g) = round(random(nogodistcountsdraw,1,1));
        end
        
        sampmotorgo = datasample(motorbounds(2):250:motorbounds(3),5000);
        sampthetago = polyval(coefftheta,sampmotorgo,[],mutheta);
        sampcountsgo = polyval(coeffcounts,sampmotorgo,[],mucounts);
        for g = 1:length(sampmotorgo)
            godistthetadraw = makedist('normal','mu',sampthetago(g),'sigma',sigvaluestheta(p));
            godistcountsdraw = makedist('normal','mu',sampcountsgo(g),'sigma',sigvaluescounts(p));
            godisttheta(g) = random(godistthetadraw,1,1);
            godistcounts(g) = round(random(godistcountsdraw,1,1));
        end
        
        valuesnogotheta = histcounts(nogodisttheta,optboundstheta(end):1:optboundstheta(1));
        valuesgotheta = histcounts(godisttheta,optboundstheta(end):1:optboundstheta(1));
        valuesnogocounts = histcounts(nogodistcounts,optboundscounts);
        valuesgocounts = histcounts(godistcounts,optboundscounts);
        
        %Plotting theta distributions
        figure(380);subplot(2,length(sigvaluestheta),p)
        bar(linspace(optboundstheta(end),optboundstheta(1),length(valuesnogotheta)),valuesnogotheta,'r');
        alpha(.5)
        hold on; bar(linspace(optboundstheta(end),optboundstheta(1),length(valuesnogotheta)),valuesgotheta,'b')
        alpha(.5)
        set(gca,'xlim',[-10 35],'ylim',[0 1500])
        
        figure(380);subplot(2,length(sigvaluestheta),p+length(sigvaluestheta))
        bar(linspace(optboundstheta(end),optboundstheta(1),length(valuesnogotheta)),valuesnogotheta,'r');
        alpha(.5)
        hold on; bar(linspace(optboundstheta(end),optboundstheta(1),length(valuesnogotheta)),valuesgotheta,'b')
        alpha(.5)
        set(gca,'xlim',[optboundstheta(2)-5 optboundstheta(2)+5],'ylim',[0 750])
        
        %Plotting count distributions
        figure(381);subplot(1,length(sigvaluestheta),p)
        bar(0:15,valuesnogocounts,'r');
        alpha(.5)
        hold on; bar(0:15,valuesgocounts,'b');
        alpha(.5)
        set(gca,'xlim',[-.5 15.5])
        
        
    end
    set(gcf, 'Units', 'pixels', 'Position', [500, 500, 2000, 400]);
    
    
    
    
    
    
    
end
end




