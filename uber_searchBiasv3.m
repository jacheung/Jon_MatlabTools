ranges = -75:1:75;
% sbiastwo = nan(length(U),2);
% sbiasthree = nan(length(U),1);
% sbiasfour = nan(length(U),1);

type = {U2};
% type = {Nx,BVx};
clearvars -except U D BV SM BVx Nx type ranges POP
close all
for q=1
    U = type{q};
    
    
    %%
    for rec =1:length(U)
        
        
        farthestnogoT= find(U{rec}.meta.motorPosition<min(U{rec}.meta.motorPosition)+10000);
        
%         %choosing trials with which to evaluate exploration
%         if strcmp(U{rec}.meta.layer,'D')
%             farthestnogoT = find(U{rec}.meta.trialType==0);
%         else
%             farthestnogoT= find(U{rec}.meta.motorPosition<min(U{rec}.meta.motorPosition)+10000);
%         end
%         
        %Only use farthest nogo to calculate thetas (uninhibited whisking)
        thetasall = squeeze(U{rec}.S_ctk(1,:,:));
        thetasnan = nan(size(thetasall));
        thetasnan(:,farthestnogoT) = 1;
        thetas = thetasall.*thetasnan;
        
%         [whisks] = findMaxMinProtraction(U{rec},'sampling');
           [whisks] = findMaxMinProtraction(U{rec},5,'sampling');
           [ALLwhisks] = findMaxMinProtraction(U{rec},2.5);
        
        selectedT=thetas(whisks.peakidx);
        selectedT(isnan(selectedT))=[];
        selectedT=sort(selectedT);
        [v,selectedT] = ecdf(selectedT);
        
        
        % FIND POLE POSITION RANGES (inherited from searchbias_v2 so may
        % have lots of excess; trim sometime... 180517)
        %p(pole at theta) distribution
        ppoles = nan(U{rec}.k,1000);
        for b = 1:U{rec}.k
            touchIdx = [find(U{rec}.S_ctk(9,:,b)==1) find(U{rec}.S_ctk(12,:,b)==1)];
            if ~isempty(touchIdx)
                tmpThetas = U{rec}.S_ctk(1,touchIdx,b);
                ppoles(b,1:length(tmpThetas)) = tmpThetas;
%                 ppoles(b,:) = [median(tmpThetas) std(tmpThetas)];
            end
        end
        
        goIdx = U{rec}.meta.trialType == 1;
        nogoIdx = U{rec}.meta.trialType == 0;
        
        gos = round([ min(ppoles(goIdx,1)) max(ppoles(goIdx,1))])     
        nogos = round([ min(ppoles(nogoIdx,1)) max(ppoles(nogoIdx,1))])

% FIND SPAN OF WHISKING       

        %Protraction and Retraction 
        matchedIdxRet = nan(length(whisks.peakidx),1);
        matchedIdxPro = nan(length(whisks.peakidx),1);
        for k = 1:length(whisks.peakidx)
            matchedtmpRet  = find(ALLwhisks.troughidx>whisks.peakidx(k),1);
            matchedtmpPro =  find(ALLwhisks.troughidx<whisks.peakidx(k),1,'last');
            if ~isempty(matchedtmpRet)
                matchedIdxRet(k) = ALLwhisks.troughidx(matchedtmpRet);
            end
            if ~isempty(matchedtmpPro)
                matchedIdxPro(k) = ALLwhisks.troughidx(matchedtmpPro);
            end
        end
        
        whiskDurationPro = whisks.peakidx - matchedIdxPro ;
        tossIdxPro = find(whiskDurationPro>60);
        spanIdxPro = [matchedIdxPro whisks.peakidx];
        spanIdxPro(tossIdxPro,:) = [];
        spanIdxPro(sum(isnan(spanIdxPro), 2) >= 1, :) = []; 
        thetaspanPro = thetas(spanIdxPro);
        spanIdxPro(sum(isnan(thetaspanPro), 2) >= 1, :) = []; 
        thetaspanPro(sum(isnan(thetaspanPro), 2) >= 1, :) = []; 
        
        whiskDurationRet = matchedIdxRet-whisks.peakidx;
        tossIdxRet = find(whiskDurationRet>60);
        spanIdxRet = [whisks.peakidx matchedIdxRet];
        spanIdxRet(tossIdxRet,:) = [];
        spanIdxRet(sum(isnan(spanIdxRet), 2) >= 1, :) = [];        
        thetaspanRet = thetas(spanIdxRet);
        spanIdxRet(sum(isnan(thetaspanRet), 2) >= 1, :) = []; 
        thetaspanRet(sum(isnan(thetaspanRet), 2) >= 1, :) = []; 
        
        
        
        
        
        
        
        
%CALCULATION PROPORTION OF GO AND NOGO RANGE CONTACTED/WHISK
        goranges = gos(1):gos(end);
        nogoranges = nogos(1):nogos(end);
        
    %Protraction 
        goprop = nan(length(thetaspanPro),1);
        nogoprop = nan(length(thetaspanPro),1);
        
        for p = 1:length(thetaspanPro)
            tmpprorange = round(thetaspanPro(p,:));
            proVals=tmpprorange(1):tmpprorange(end);
            gointersect = intersect(goranges,proVals);
            nogointersect = intersect(nogoranges,proVals);
            
            if ~isempty(gointersect)
                goprop(p) = numel(gointersect)./numel(goranges);
            else
                goprop(p) = 0;
            end
            
            if ~isempty(gointersect)
                nogoprop(p) = numel(nogointersect)./numel(nogoranges);
            else
                nogoprop(p) = 0;
            end
        end
        
        protractionEBias = (nogoprop - goprop)./(nogoprop + goprop);
     %Retraction 
        gopropRet = nan(length(thetaspanRet),1);
        nogopropRet = nan(length(thetaspanRet),1);
        
         for p = 1:length(thetaspanRet)
            tmpretrange = round(thetaspanRet(p,:));
            retVals=tmpretrange(end):tmpretrange(1);
            gointersect = intersect(goranges,retVals);
            nogointersect = intersect(nogoranges,retVals);
            
            if ~isempty(gointersect)
                gopropRet(p) = numel(gointersect)./numel(goranges);
            else
                gopropRet(p) = 0;
            end
            
            if ~isempty(gointersect)
                nogopropRet(p) = numel(nogointersect)./numel(nogoranges);
            else
                nogopropRet(p) = 0;
            end
         end
        
        retractionEBias = (nogopropRet - gopropRet)./(nogopropRet + gopropRet); 
        
        
        totalBias{q}{rec} = [protractionEBias;retractionEBias]; 
        U{rec}.meta.sbias = nanmean([protractionEBias;retractionEBias]);
        
%       
%find vals for a single trial 
%          protractionEBias(find(ceil(spanIdxPro(:,1)/4000)==148))
%          retractionEBias(find(ceil(spanIdxRet(:,1)/4000)==148))
      
%         totalBias{q}{rec} = [protractionEBias]; 
%         totalBias{q}{rec} = [retractionEBias]; 

    end
end  
%% Some plotting stuff
colors = {'DarkGreen','DarkMagenta','DarkTurquoise'};
semivals = cellfun(@nanmean,totalBias{1});
contvals = cellfun(@nanmean,totalBias{2});

POP.SBIAS{2} = semivals;
POP.SBIAS{3} = contvals;

figure(480);clf;
h=scatter(semivals,ones(length(semivals),1).*2,'filled')
hold on; errorbar(nanmean(semivals),2,nanstd(semivals),'horizontal','k');
 h.CData = rgb(colors{2});
h = scatter(contvals,ones(length(semivals),1).*1,'filled');
 h.CData = rgb(colors{3});
hold on;errorbar(nanmean(contvals),1,nanstd(contvals),'horizontal','k');
set(gca,'xlim',[-1 1],'xtick',[-1:1:1],'ylim',[.5 2.5],'ytick',[]) 
xlabel('Exploration bias');

%plotting individual mouse
selectedBias = cell2mat(totalBias{1}');
sortedvals = binslin(selectedBias,selectedBias,'equalE',22,-1,1);
counts = cellfun(@numel,sortedvals);
proportions = counts./sum(counts);
figure(3280);
hold on; bar(-1:.1:1,proportions,'facecolor',rgb(colors{2}))
set(gca,'xlim',[-1.05 1.05],'xtick',[-1:1:1],'ytick',[0:.1:1],'ylim',[0 1])

% ANOVA FOR SIG
alpha = .01; %significance level of .01 
[p,tbl,stats] = anova1([semivals' contvals']); %anova1 
comp = multcompare(stats); %comparison between all groups. 
bonfcorr = alpha/length(contvals); %post hoc bonferroni correction of pval alpha/n

[comp(:,1:2) comp(:,end)<bonfcorr ]
        
        
 