type = {SM,BV};


collat =[];thetacollat = [];

for d = 2
    colors = {'DarkMagenta','DarkTurquoise'};
    U=BV
    
    gocountssum = zeros(1,16);
    nogocountssum = zeros(1,16);
    
    clear g;clear ng
    [V] = classifierWrapper(U);
    for rec = 1:length(U)
        nogotouches = [V(rec).var.FA{1} V(rec).var.CR{1}]; %makes sure there is at least one touch before using array for building distribution
        
        if ~isempty(nogotouches) || numel(nogotouches)>1
            gos = [V(rec).touchNum.hit V(rec).touchNum.miss];
            nogos = [V(rec).touchNum.CR V(rec).touchNum.FA];
            totalT=numel([gos nogos]);
            goDist = histc(gos,0:15);
            nogoDist = histc(nogos,0:15);
            
            %LICK PROBABILITY
            golick = [V(rec).touchNum.hit];
            gonolick =  [V(rec).touchNum.miss];
            nogolick = [V(rec).touchNum.FA];
            nogonolick = [V(rec).touchNum.CR];
            
            [glcounts] = histc(golick,0:15);
            [gnlcounts] = histc(gonolick,0:15);
            [nlcounts] = histc(nogolick,0:15);
            [ngnlcounts] = histc(nogonolick,0:15);
            
            golickp = glcounts./ (glcounts+gnlcounts);
            nogolickp = nlcounts./ (nlcounts+ngnlcounts);
            grouplickp = (glcounts+nlcounts) ./ (glcounts+gnlcounts+nlcounts+ngnlcounts);
            
            gocountssum = gocountssum+glcounts;
            nogocountssum = nogocountssum + nlcounts;
 
            figure(50);subplot(2,5,rec);
            bar(0:15,goDist./numel(gos),'b')
            alpha(.5)
            hold on; bar(0:15,nogoDist./numel(nogos),'r');
            alpha(.5)
            set(gca,'xlim',[-.5 16],'xtick',0:5:15,'ylim',[0 1],'ytick',[0 .5 1])
            hold on; plot(0:15,golickp,'b')
            hold on; plot(0:15,nogolickp,'r')
            
            glick{rec}=golickp;
            nglick{rec}=nogolickp;
            grouplick{rec} = grouplickp;
            
            g{rec}=goDist./numel(gos);
            ng{rec}=nogoDist./numel(nogos);
            
            gngratio(rec) = sum(goDist.*(0:15))./sum(nogoDist.*(0:15));
            
        end
    end
    
    %Plotting mean gonogo distributions
    GgoDist = mean(cell2mat(g'));
    GngDist = mean(cell2mat(ng'));
    gotouch = [];
    nogotouch = [];
    
    figure(50+d);
    bar(0:15,GgoDist,'b');
    hold on; bar(0:15,GngDist,'r');
    set(gca,'xlim',[-.5 16],'xtick',0:5:15,'ylim',[0 1.1*.75],'ytick',[0 .25 .5 .75])
    set(gcf, 'Units', 'pixels', 'Position', [0, 0, 500, 600]);
    %plotting lick probabilty ontop of means
    golmeans = nanmean(cell2mat(glick'));
    nogolmeans = nanmean(cell2mat(nglick'));
    gocibin=zeros(2,length(gocountssum));
    nogocibin=zeros(2,length(nogocountssum));%doing 95% ci of each bin in relation to each bin
    for i=1:length(gocountssum)
%         x=groupKappgo(:,i)
%         SEM = nanstd(x)/sqrt(length(x));               % Standard Error
%         ts = tinv([0.025  0.975],length(x)-1);      % T-Score
%         Tcibingo([2 1],i) = nanmean(x) + ts*SEM;   
%        
        [~, pci] = binofit(golmeans(i).*gocountssum(i),gocountssum(i));
        gocibin([2 1],i) = pci;
        [~, pcing] = binofit(nogolmeans(i).*nogocountssum(i),nogocountssum(i));
        nogocibin([2 1],i) = pcing;
    end
    
    hold on; plot(0:15,nanmean(cell2mat(glick'))*.75,'b','linewidth',3)
    hold on; plot(0:15,gocibin(1,:)*.75,'b');
    hold on; plot(0:15,gocibin(2,:)*.75,'b');
    hold on; plot(0:15,nanmean(cell2mat(nglick'))*.75,'r','linewidth',3)
    hold on; plot(0:15,nogocibin(1,:)*.75,'r');
    hold on; plot(0:15,nogocibin(2,:)*.75,'r');
    
    %INDIVIDUAL OVERLAP BETWEEEN GO AND NOGO;
    INDIVgoproportionsOG = cell2mat(g')./sum(cell2mat(g'),2);
    INDIVnogoproportionsOG = cell2mat(ng')./sum(cell2mat(ng'),2);
    overlapIDX = double(INDIVgoproportionsOG>0).*double(INDIVnogoproportionsOG>0);
    sampledBins = double((double(INDIVgoproportionsOG>0)+double(INDIVnogoproportionsOG>0))>=1);
    
    overlapspace = sum(min(INDIVnogoproportionsOG,INDIVgoproportionsOG).*overlapIDX,2);
    totaluniquespace = sum(max(INDIVnogoproportionsOG,INDIVgoproportionsOG),2);
    countsoverlap = overlapspace./totaluniquespace;
    
    
    % CALCULATING DIFFERENCE IN LICK PROB IN OVERLAPPING BINS
    allgolicks = cell2mat(glick').*overlapIDX;
    allnglicks = cell2mat(nglick').*overlapIDX;
    meanlickDiff = nanmean(allgolicks-allnglicks,2);
    
    
    %Plotting go touches : nogotouches ratio
    allgotouches = sum(cell2mat(g').*repmat(0:15,length(g),1),2);
    allnogotouches = sum(cell2mat(ng').*repmat(0:15,length(ng),1),2);
    ratio = ((allgotouches - allnogotouches) ./(allgotouches+allnogotouches))';
    figure(430);hold on;
    h = scatter(ratio,ones(length(ratio),1)*d,'filled');
    h.CData = rgb(colors{d});
    errorbar(mean(ratio),d,std(ratio),'horizontal','ko','markersize',15,'markerfacecolor',rgb(colors{d}))
    plot([0 0],[.5 3.5],'-.k')
    set(gca,'ydir','reverse','xlim',[-1 1],'xtick',[-1 0 1],'ytick',[])
    xlabel('Go touch count to nogo touch count')
    
    
    
    gnglickProb{d} = [nanmean(cell2mat(glick')) ;nanmean(cell2mat(nglick'))];
    gngcountsratios{d} = ratio;
    POPcountsoverlap{d} = countsoverlap;
    
    sum(goDist.*(0:15))./sum(nogoDist.*(0:15))
    
    
    %% THETA DIST for one single mouse
    
    
    range = [-25:50]
    gos = [];
    nogos = [];
    groupgorange = zeros(length(U),2);
    groupnogorange = zeros(length(U),2);
    popvalsgo = zeros(length(U),length(range));
    popvalsnogo = zeros(length(U),length(range));
    for var = 1
        if var == 5
            range = linspace(-3.14,3.14,15);
        else
            range = [-25:50];
        end
        for rec = 1:length(U)
            tmpgos = [V(rec).var.hit{var} V(rec).var.miss{var}];
            tmpnogos = [V(rec).var.FA{var} V(rec).var.CR{var}];
            
            tmpgothetas = histc([tmpgos],range);
            tmpnogothetas = histc([tmpnogos],range);
            figure(32); subplot(2,5,rec);
            plot(range,tmpgothetas,'b');
            hold on; plot(range,tmpnogothetas,'r');
            
            thetagovals = range(tmpgothetas>0);
            thetanogovals = range(tmpnogothetas>0);
            groupgorange(rec,1:2) = [min(thetagovals) max(thetagovals)];
            groupnogorange(rec,1:2) = [min(thetanogovals) max(thetanogovals)];
            
            popvalsgo(rec,:) = tmpgothetas;
            popvalsnogo(rec,:) =  tmpnogothetas;
            
            
            
            if var == 5
                set(gca,'xlim',[min(range) max(range)],'xtick',linspace(-3.14,3.14,3),'xticklabel',{'-\pi','0','\pi'})
            end
            
            gos = [gos V(rec).var.hit{var} V(rec).var.miss{var}];
            nogos = [nogos V(rec).var.CR{var} V(rec).var.FA{var}];
            
            gthetas{rec}=tmpgos;
            ngthetas{rec}=tmpnogos;
            
            if rec == 1 || rec ==6
                ylabel('Number of Touches')
            elseif rec == 3 || rec == 8
                xlabel([V(1).varNames{var} ' at Touch'])
            end
            
            
        end
        
        gothetas = histc([gos],range);
        nogothetas = histc([nogos],range);
        figure(89+d); clf;
        bar(range,gothetas./(sum(gothetas)+sum(nogothetas)),'b');
        hold on;bar(range,nogothetas./(sum(gothetas)+sum(nogothetas)),'r');
        if var == 5
            hold on; bar(range,nogothetas./sum(nogothetas),'r');
            set(gca,'xlim',[min(range) max(range)],'xtick',linspace(-3.14,3.14,5),'xticklabel',{'-\pi','-\pi/2','0','\pi/2','\pi'})
        end
        set(gca,'ylim',[0 .05],'ytick',0:.025:.05,'xtick',-40:20:60,'xlim',[-35 50])
%         ylabel('Proportion of touches')
%         xlabel([V(1).varNames{var} ' at touch'])
        set(gcf, 'Units', 'pixels', 'Position', [0, 0, 500, 600]);
        
        
    end
    
    
    %THETA OVERLAP
    INDIVgoproportions = popvalsgo./sum(popvalsgo,2);
    INDIVnogoproportions = popvalsnogo./sum(popvalsnogo,2);
    overlapIDX = double(INDIVgoproportions>0).*double(INDIVnogoproportions>0);
    INDIVgoproportions = INDIVgoproportions.*overlapIDX;
    INDIVnogoproportions = INDIVnogoproportions.*overlapIDX;
    
    
    difference = INDIVnogoproportions - INDIVgoproportions;
    overlap1keep = INDIVgoproportions.*double(difference>0);
    pt2keep = INDIVnogoproportions .* double(difference<=0);
    overlaptmp = sum([overlap1keep pt2keep],2);
    thetaoverlap = overlaptmp./2;
    
    POPthetaoverlap{d} = thetaoverlap;
    
    
    %Plotting median of theta values +/- 2*standard deviation for full range of
    %theta values
    figure(809); hold on
    errorbar(cellfun(@median,ngthetas),(length(ngthetas)*d-length(ngthetas)+1:length(ngthetas)*d),2*cellfun(@std,ngthetas),'horizontal','ro','markersize',10,'markerfacecolor',rgb(colors{d}))
    errorbar(cellfun(@median,gthetas),(length(gthetas)*d-length(gthetas)+1:length(gthetas)*d),2*cellfun(@std,gthetas),'horizontal','bo','markersize',10,'markerfacecolor',rgb(colors{d}))
    set(gca,'ytick',[],'ylim',[0 31],'ydir','reverse')
    set(gcf, 'Units', 'pixels', 'Position', [0, 0, 500, 600]);
    %calculating task difficulty: ( median(ng) - 2std(ng)) - (median(g)+2std(g))
    gngtaskease{d} = (cellfun(@median,ngthetas)-(2*cellfun(@std,ngthetas)))   - (cellfun(@median,gthetas)+(2*cellfun(@std,gthetas)));
    %calculating means of theta ranges
    gngfullthetaranges{d} = (cellfun(@max,ngthetas)-cellfun(@min,ngthetas)) + (cellfun(@max,gthetas)-cellfun(@min,gthetas));
    
    gowidth = groupgorange(:,2)-groupgorange(:,1);
    nogowidth = groupnogorange(:,2)-groupnogorange(:,1);
    
    
    % figure(321);hold on;h=scatter(iv{d},thetameandiff,'filled');
    % h.CData = rgb(colors{d});
    
    
end

%% Plotting distribution of features
colors = {'DarkMagenta','DarkTurquoise'};
xdata = POPthetaoverlap;
% xdata = POPcountsoverlap;

for d = 1:length(xdata)
    figure(680);hold on
    h = scatter(xdata{d},ones(length(xdata{d}),1).*d,'filled');
    h.CData = rgb(colors{d});
    errorbar(mean(xdata{d}),d,std(xdata{d}),'horizontal','ko','markersize',15,'markerfacecolor',rgb(colors{d}))
    
end

set(gca,'ydir','reverse','ylim',[.5 2.5],'ytick',[],'xlim',[0 .5],'xtick',[0:.25:1])
set(gcf, 'Units', 'pixels', 'Position', [0, 0, 400, 600]);
% ANOVA FOR SIG
alpha = .01; %significance level of .01
[p,tbl,stats] = anova1(cell2mat(xdata)); %anova1
comp = multcompare(stats);

%% Plotting for SBIASxtouch count distribution figure320
% Plotting for SBIAS x meanof theta distributoin figure 321
xdata = [POPv2.SBIAS{2:3}];
% ydata = cell2mat(gngcountsratios);
ydata = cell2mat(POPcountsoverlap')';
thetaydata = cell2mat(gngfullthetaranges);
[~, orders] = sort(xdata);

[coeff, ~ , mu] = polyfit(xdata(orders),ydata(orders),1);
f = polyval(coeff,xdata(orders),[],mu);

[tcoeff, ~ , tmu] = polyfit(xdata(orders),thetaydata(orders),1);
t = polyval(tcoeff,xdata(orders),[],tmu);


modelvals = fitlm(xdata,ydata);
tmodelvals = fitlm(xdata,thetaydata);

pval = modelvals.Coefficients{2,4};
adjrsq = modelvals.Rsquared.Adjusted;

tpval = tmodelvals.Coefficients{2,4};
tadjrsq = tmodelvals.Rsquared.Adjusted;

figure(320); plot(xdata(orders),f,'k');
legend('Semi-Continuous','Continuous','location','northwest')
set(gca,'ylim',[-1 1],'ytick',[-1:.5:1],'xlim',[-5 25])
disp(['Rsquared = ' num2str(adjrsq) ' and pvalue = ' num2str(pval)])
xlabel('Search Bias')
ylabel('More nogo touches ------ More go touches')
set(figure(320), 'Units', 'pixels', 'Position', [0, 0, 600, 750]);

figure(321); plot(xdata(orders),t,'k');
legend('Semi-Continuous','Continuous','location','northwest')
set(gca,'xlim',[-5 30])
disp(['Rsquared = ' num2str(tadjrsq) ' and pvalue = ' num2str(tpval)])
xlabel('Gap (theta)')
ylabel('Difference between mean of go - mean of nogo')
set(figure(321), 'Units', 'pixels', 'Position', [0, 0, 600, 750]);


%% Task D x SBIAS CORRELATIOn


xdata = [POPv2.SBIAS{2:3}]';
% xdata = cell2mat(gngtaskease);

% ydata = cell2mat(gngcountsratios)*-1
% ydata = cell2mat(gngfullthetaranges);
ydata = cell2mat(POPthetaoverlap');
[~, orders] = sort(xdata);

[coeff, ~ , mu] = polyfit(xdata(orders),ydata(orders),1);
g = polyval(coeff,xdata(orders),[],mu);

modelvals = fitlm(xdata,ydata);
pval = modelvals.Coefficients{2,4};
adjrsq = modelvals.Rsquared.Adjusted;
slope = modelvals.Coefficients{2,1};

disp(['Rsquared = ' num2str(adjrsq) ' and pvalue = ' num2str(pval) ' and slope = ' num2str(slope)])



colors = {'DarkMagenta','DarkTurquoise'};
ranges = [1:10;11:20];
figure(23);clf;
for d = 1:length(colors)
    hold on; h=scatter(xdata(ranges(d,:)),ydata(ranges(d,:)),'filled');
    h.CData = rgb(colors{d});
end
% hold on; plot([-1 1],[0 0],'-.k')
hold on; plot(xdata(orders),g,'k')

set(gca,'ylim',[0 .5],'ytick',-1:.5:1,'xlim',[-1 0],'xtick',[-1:1:1])

set(figure(23), 'Units', 'pixels', 'Position', [0, 0, 600, 600]);
legend('Semi-Continuous','Continuous','location','northwest')


