type = {SM,BV};

for d = 1
    U=type{d};
 countgroupvals = [];
 polegroupvals = [];
 groupcountpoles = cell(16,1);
    for rec = 1:length(U)
        [prelickTouches] = assist_predecisionVar(U{rec});
        cpoles = cellfun(@numel,prelickTouches)';
        poles = U{rec}.meta.motorPosition;
        poleranges = U{rec}.meta.ranges;
        gos = find(U{rec}.meta.trialType==1);
        nogos = find(U{rec}.meta.trialType==0);

        
        %LICKS
        hits = double(U{rec}.meta.trialCorrect == 1) .* double(U{rec}.meta.trialType ==1);
        miss = double(U{rec}.meta.trialCorrect == 0) .* double(U{rec}.meta.trialType ==1);
        FA = double(U{rec}.meta.trialCorrect == 0) .* double(U{rec}.meta.trialType ==0);
        CR = double(U{rec}.meta.trialCorrect == 1) .* double(U{rec}.meta.trialType ==0);
        lix = hits + FA;
        correct = U{rec}.meta.trialCorrect;

        
        against = lix; 
        poles = poles;
        cpoles = cpoles; 
        
        %COUNTS x POLE POSITIONS 
        normpoles = normalize_var(poles,-1,1);
        countspolessorted = binslin(cpoles,normpoles','equalE',17,-.5,15.5);
        meanpoles = cellfun(@nanmean,countspolessorted);
        stdpoles = cellfun(@nanstd,countspolessorted);
        figure(4700);subplot(2,5,rec)
        errorbar(0:15,meanpoles,stdpoles,'-ko')
%         hold on; plot([-1 16],[mean(poleranges) mean(poleranges)],'-.k')
        hold on; plot([-1 16],[0 0],'-.k')
         set(gca,'Ydir','reverse','ylim',[-1 1],'ytick',[],'xtick',0:5:15,'xlim',[-1 16])
        
        for g = 1:length(countspolessorted)
            groupcountpoles{g} = [groupcountpoles{g};countspolessorted{g}];
        end
        
         keepIDX = intersect(find(cpoles>=2),find(cpoles<=4));
        
        
        %Bin based on pole positions or number of touches
        [sortedCounts]= binslin(poles,[cpoles against' U{rec}.meta.trialType'],'equalE',11,poleranges(1),poleranges(end));
        [sortedPoles]= binslin(cpoles,[poles' against' U{rec}.meta.trialType'],'equalE',17,-.5,15.5);
        
        %Find means of each bin to mean normalize
        polesmeans = cell2mat(cellfun(@(x) nanmean(x(:,1)),sortedPoles,'uniformoutput',0));
        countsmeans = cell2mat(cellfun(@(x) nanmean(x(:,1)),sortedCounts,'uniformoutput',0));

        
        allcountvals = [];
        allpolevals = [];
        for p = 1:length(sortedCounts) 
            gvals = [sortedCounts{p}(:,1)-countsmeans(p) sortedCounts{p}(:,2:3)];
            allcountvals = [allcountvals; gvals];
        end
        
        for j = 1:length(sortedPoles) 
            gvals = [sortedPoles{j}(:,1)-polesmeans(j) sortedPoles{j}(:,2:3)];
            gvalstmp{j} = gvals;
            allpolevals = [allpolevals;gvals];
        end
        
        %CROSSING IDX
        crossbounds = [meanpoles-(.5*stdpoles) meanpoles+(.5*stdpoles)];
        crossidx = intersect(find(crossbounds(:,1)<=0),find(crossbounds(:,2)>=0));
        
        for t = 1:length(gvalstmp)
            allpolevalstmp = gvalstmp{t};
            [sortedpoleVals] =binslin(allpolevalstmp(:,1),allpolevalstmp(:,2),'equalE',12,-50000,50000);
            lickPs = cell2mat(cellfun(@(x) nanmean(x,1),sortedpoleVals,'uniformoutput',0));
            figure(4800);hold on; subplot(2,5,rec)
            
            if ~isempty(intersect(t,crossidx))
                plot(linspace(-1,1,length(lickPs))+meanpoles(t),lickPs+t,'-ro')
            else
            plot(linspace(-1,1,length(lickPs))+meanpoles(t),lickPs+t,'-ko')
            end
            set(gca,'xlim',[-1.1 1.1],'xtick',-1:1:1,'ylim',[0 10],'ytick',0:5:10)
        end
        
        
        % Plotting for individual mice holding POLES
        [sortedVals ] =binslin(allcountvals(:,1),allcountvals(:,2),'equalE',22,-10,10);
        lickPs = cell2mat(cellfun(@mean,sortedVals,'uniformoutput',0));
        figure(500);subplot(2,5,rec); hold on;
        plot(linspace(-10,10,length(lickPs)),lickPs,'-ko')
        set(gca,'xlim',[-10 10],'xtick',-10:10:10,'ylim',[0 1.1],'ytick',0:.5:1)
        % Plotting for individual mice holding poles GONOGO
        govals = allcountvals(:,3)==1;
        nogovals = allcountvals(:,3)==0;
        [sortedVals ] =binslin(allcountvals(govals,1),allcountvals(govals,2),'equalE',22,-10,10);
        lickPs = cell2mat(cellfun(@mean,sortedVals,'uniformoutput',0));
        figure(501);subplot(2,5,rec); hold on;
        plot(linspace(-10,10,length(lickPs)),lickPs,'-bo')
        set(gca,'xlim',[-10 10],'xtick',-10:10:10,'ylim',[0 1.1],'ytick',0:.5:1)
        [sortedVals ] =binslin(allcountvals(nogovals,1),allcountvals(nogovals,2),'equalE',22,-10,10);
        lickPs = cell2mat(cellfun(@mean,sortedVals,'uniformoutput',0));
        figure(501);subplot(2,5,rec); hold on;
        plot(linspace(-10,10,length(lickPs)),lickPs,'-ro')
        set(gca,'xlim',[-10 10],'xtick',-10:10:10,'ylim',[0 1.1],'ytick',0:.5:1)

        
%         % Plotting for individual mice holding COUNTS
%         [sortedpoleVals] =binslin(allpolevals(:,1),allpolevals(:,2),'equalE',12,-50000,50000);
%         lickPs = cell2mat(cellfun(@(x) nanmean(x,1),sortedpoleVals,'uniformoutput',0));
%         figure(380);subplot(2,5,rec); hold on;
%         plot(linspace(-5,5,length(lickPs)),lickPs,'-ko')
%         set(gca,'xlim',[-5 5],'xtick',-5:5:5,'ylim',[0 1.1],'ytick',0:.5:1)
%         % Plotting for individual mice holding COUNTS GONOGO
%         govals = allpolevals(:,3)==1;
%         nogovals = allpolevals(:,3)==0;
%         [sortedpoleVals] =binslin(allpolevals(govals,1),allpolevals(govals,2),'equalE',12,-50000,50000);
%         lickPs = cell2mat(cellfun(@(x) nanmean(x,1),sortedpoleVals,'uniformoutput',0));
%         figure(381);subplot(2,5,rec); hold on;
%         plot(linspace(-5,5,length(lickPs)),lickPs,'-bo')
%         set(gca,'xlim',[-5 5],'xtick',-5:5:5,'ylim',[0 1.1],'ytick',0:.5:1)
%          [sortedpoleVals] =binslin(allpolevals(nogovals,1),allpolevals(nogovals,2),'equalE',12,-50000,50000);
%         lickPs = cell2mat(cellfun(@(x) nanmean(x,1),sortedpoleVals,'uniformoutput',0));
%         figure(381);subplot(2,5,rec); hold on;
%         plot(linspace(-5,5,length(lickPs)),lickPs,'-ro')
%         set(gca,'xlim',[-5 5],'xtick',-5:5:5,'ylim',[0 1.1],'ytick',0:.5:1)
        
        %keeping raw values for group analysis
        countgroupvals = [countgroupvals ;allcountvals];
        polegroupvals = [polegroupvals ;allpolevals];
           
    end
    
    
    %% POPULATION PLOTTING COUNTS x POSITION
    
    yax = cellfun(@mean,groupcountpoles) ;
    yaxstd = cellfun(@std,groupcountpoles);
    figure(4888);clf;errorbar(0:15,yax,yaxstd,'-ko')
    hold on;plot([-1 16],[0 0],'-.k')
    set(gca,'Ydir','reverse','ylim',[-1 1],'ytick',[],'xtick',0:5:15)
    
    
    %% POPULATION PLOTTING  POLES
%         [sortedpoleVals] =binslin(polegroupvals(:,1),polegroupvals(:,2),'equalE',12,-50000,50000);
        [sortedpoleVals,sortedBypoleVals] =binslin(polegroupvals(:,1),polegroupvals(:,2),'equalN',10);
        lickPs = cell2mat(cellfun(@(x) nanmean(x(:,1)),sortedpoleVals,'uniformoutput',0));
        xax = cellfun(@mean,sortedBypoleVals);
        Tcibin=zeros(2,length(sortedpoleVals)); %doing 95% ci of each bin in relation to each bin
        for i=1:length(sortedpoleVals)
        [phat, pci] = binofit(sum(sortedpoleVals{i}),numel(sortedpoleVals{i}));
        Tcibin([2 1],i) = pci;
        end    
        figure(383);clf
        plot(xax./10000,lickPs,'-ko')
        hold on;plot(xax./10000,Tcibin(1,:),'-.k')
        hold on;plot(xax./10000,Tcibin(2,:),'-.k')
        set(gca,'xlim',[-5 5],'xtick',-5:5:5,'ylim',[0 1.1],'ytick',0:.5:1)
        
        govals = polegroupvals(:,3)==1;
        nogovals = polegroupvals(:,3)==0;
        
%       [sortedpoleVals] =binslin(polegroupvals(govals,1),polegroupvals(govals,2),'equalE',12,-50000,50000);
        [sortedpoleValsgo,sortedByGo] =binslin(polegroupvals(govals,1),polegroupvals(govals,2),'equalN',10);
        lickPs = cell2mat(cellfun(@(x) nanmean(x,1),sortedpoleValsgo,'uniformoutput',0));
        xaxgo = cellfun(@mean,sortedByGo);
        Tcibin=zeros(2,length(sortedpoleValsgo)); %doing 95% ci of each bin in relation to each bin
        for i=1:length(sortedpoleValsgo)
        [phat, pci] = binofit(sum(sortedpoleValsgo{i}),numel(sortedpoleValsgo{i}));
        Tcibin([2 1],i) = pci;
        end    
        figure(384);clf
        plot(xaxgo./10000,lickPs,'-bo')
        hold on;plot(xaxgo./10000,Tcibin(1,:),'-.b')
        hold on;plot(xaxgo./10000,Tcibin(2,:),'-.b')
        
%       [sortedpoleVals] =binslin(polegroupvals(nogovals,1),polegroupvals(nogovals,2),'equalE',12,-50000,50000);
        [sortedpoleValsnogo,sortedByNoGo] =binslin(polegroupvals(nogovals,1),polegroupvals(nogovals,2),'equalN',10);       
        lickPs = cell2mat(cellfun(@(x) nanmean(x,1),sortedpoleValsnogo,'uniformoutput',0));
        xaxnogo = cellfun(@mean,sortedByNoGo);
        Tcibin=zeros(2,length(sortedpoleValsnogo)); %doing 95% ci of each bin in relation to each bin
        for i=1:length(sortedpoleValsnogo)
        [phat, pci] = binofit(sum(sortedpoleValsnogo{i}),numel(sortedpoleValsnogo{i}));
        Tcibin([2 1],i) = pci;
        end    
        
        hold on; plot(xaxnogo./10000,lickPs,'-ro')
        hold on;plot(xaxnogo./10000,Tcibin(1,:),'-.r')
        hold on;plot(xaxnogo./10000,Tcibin(2,:),'-.r')
        set(gca,'xlim',[-5 5],'xtick',-5:5:5,'ylim',[0 1.1],'ytick',0:.5:1)
        
      
    %POPULATION PLOTTING COUNTS       
        [sortedVals] =binslin(countgroupvals(:,1),countgroupvals(:,2),'equalE',22,-10,10);
        lickPs = cell2mat(cellfun(@mean,sortedVals,'uniformoutput',0));
        Tcibin=zeros(2,length(sortedVals)); %doing 95% ci of each bin in relation to each bin
        for i=1:length(sortedVals)
        [phat, pci] = binofit(sum(sortedVals{i}),numel(sortedVals{i}));
        Tcibin([2 1],i) = pci;
        end
        
        figure(502);
        plot(linspace(-10,10,length(lickPs)),lickPs,'-ko')
        hold on;plot(linspace(-10,10,length(lickPs)),Tcibin(1,:),'-.k')
        hold on;plot(linspace(-10,10,length(lickPs)),Tcibin(2,:),'-.k')
        xlabel('counts') 
         set(gca,'xlim',[-10 10],'xtick',-10:5:10,'ylim',[0 1.1],'ytick',0:.5:1)

        
%         figure;boundedline(linspace(-10,10,length(lickPs)),repmat(lickPs,1,2)',Tcibin')
%         set(gca,'xlim',[-10 10],'xtick',-10:10:10,'ylim',[0 1.1],'ytick',0:.5:1)

        govals =countgroupvals(:,3)==1;
        nogovals = countgroupvals(:,3)==0;
        
        %GOS
        [sortedVals ] =binslin(countgroupvals(govals,1),countgroupvals(govals,2),'equalE',22,-10,10);
        lickPs = cell2mat(cellfun(@mean,sortedVals,'uniformoutput',0));      
        Tcibin=zeros(2,length(sortedVals)); %doing 95% ci of each bin in relation to each bin
        for i=1:length(sortedVals)
        [phat, pci] = binofit(sum(sortedVals{i}),numel(sortedVals{i}));
        Tcibin([2 1],i) = pci;
        end
        figure(504);
        plot(linspace(-10,10,length(lickPs)),lickPs,'-b','linewidth',3)
        hold on;plot(linspace(-10,10,length(lickPs)),Tcibin(1,:),'-b')
        hold on;plot(linspace(-10,10,length(lickPs)),Tcibin(2,:),'-b')
     xlabel('counts') 
        %NOGOS
        [sortedVals ] =binslin(countgroupvals(nogovals,1),countgroupvals(nogovals,2),'equalE',22,-10,10);
        lickPs = cell2mat(cellfun(@mean,sortedVals,'uniformoutput',0));
        Tcibin=zeros(2,length(sortedVals)); %doing 95% ci of each bin in relation to each bin
        for i=1:length(sortedVals)
        [phat, pci] = binofit(sum(sortedVals{i}),numel(sortedVals{i}));
        Tcibin([2 1],i) = pci;
        end
        hold on;plot(linspace(-10,10,length(lickPs)),lickPs,'-r','linewidth',3)
        hold on;plot(linspace(-10,10,length(lickPs)),Tcibin(1,:),'-r')
        hold on;plot(linspace(-10,10,length(lickPs)),Tcibin(2,:),'-r')
        set(gca,'xlim',[-5 5],'xtick',-10:5:10,'ylim',[0 1.1],'ytick',0:.5:1)
    
end
            
            
        
 