type = {D,SM,BV};

% type = {BV};
figure(320);clf;
figure(321);clf;
collat =[];thetacollat = [];
iv = POP.SBIAS;
for d = 1:length(type) 
colors = {'DarkGreen','DarkMagenta','DarkTurquoise'};
    U=type{d} ;

clear g;clear ng
[V] = classifierWrapper(U);
for rec = 1:length(U)
    nogotouches = [V(rec).var.FA{1} V(rec).var.CR{1}]; %makes sure there is at least one touch before using array for building distribution
    
    if ~isempty(nogotouches) || numel(nogotouches)>1
        gos = [V(rec).touchNum.hit V(rec).touchNum.miss];
        nogos = [V(rec).touchNum.CR V(rec).touchNum.FA];
        totalT=numel([gos nogos]);
        goDist = histc(gos,0:10);
        nogoDist = histc(nogos,0:10);
        
        figure(50);subplot(2,5,rec);
        bar(0:10,goDist./numel(gos),'b');alpha(.35)
        hold on; bar(0:10,nogoDist./numel(nogos),'r');alpha(.35)
        set(gca,'xlim',[-.5 11],'xtick',0:5:10,'ylim',[0 .5],'ytick',[0 .25 .5])
        
        g{rec}=goDist;
        ng{rec}=nogoDist;
        
        gtrimmed{rec}=goDist(2:end);
        ngtrimmed{rec}=nogoDist(2:end);
    end
end

GgoDist = sum(cell2mat(g'));
GngDist = sum(cell2mat(ng'));
gotouch = [];
nogotouch = [];
for i = 0:10
    gotouch = [gotouch ; repmat(i,GgoDist(i+1),1)];
    nogotouch = [nogotouch; repmat(i,GngDist(i+1),1)];
end

gos = [mean(gotouch) std(gotouch)];
nogos = [mean(nogotouch) std(nogotouch)];

figure(50+d);
bar(0:10,GgoDist./sum(GgoDist),'b');alpha(.35);
hold on; bar(0:10,GngDist./sum(GngDist),'r');alpha(.35);
set(gca,'xlim',[-.5 11],'xtick',0:5:10,'ylim',[0 .5],'ytick',[0 .25 .5])
% print(figure(51),'-dtiff',['Z:\Users\Jon\Projects\Characterization\' U{rec}.meta.layer '\Figures\' U{rec}.meta.layer '_countsDistribution' ])
% print(figure(50),'-dtiff',['Z:\Users\Jon\Projects\Characterization\' U{rec}.meta.layer '\Figures\' U{rec}.meta.layer '_countsDistributionINDIV' ])
%

gngrats=(cellfun(@sum,gtrimmed)-cellfun(@sum,ngtrimmed))./(cellfun(@sum,ngtrimmed)+cellfun(@sum,gtrimmed));
figure(320);hold on;h=scatter(iv{d},gngrats,'filled');
h.CData = rgb(colors{d});
set(gca,'ylim',[-1 1],'ytick',[-1:.5:1])

collat = [collat gngrats];

%% THETA DIST for one single mouse


% gos = [V(rec).var.hit{1} V(rec).var.miss{1}];
% nogos = [V(rec).var.CR{1} V(rec).var.FA{1}];
% thetas = histc([gos nogos],[-50:50]);
% figure(23);subplot(2,1,1)
% plot([-50:50],thetas./numel([gos nogos]),'color',[.5 .5 .5]);
% set(gca,'xtick',[-50:25:50],'xticklabel',[-50:25:50]);
% title('Theta Distribution')
% ylabel('Proportion of Touches')


gos = [];
nogos = [];
groupgorange = zeros(length(U),2);
groupnogorange = zeros(length(U),2);
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
        plot(range,tmpgothetas,'b');alpha(.5)
        hold on; plot(range,tmpnogothetas,'r');alpha(.5)
        
        thetagovals = range(tmpgothetas>0);
        thetanogovals = range(tmpnogothetas>0);
        groupgorange(rec,1:2) = [min(thetagovals) max(thetagovals)];
        groupnogorange(rec,1:2) = [min(thetanogovals) max(thetanogovals)];
        
        
        
        if var == 5
            set(gca,'xlim',[min(range) max(range)],'xtick',linspace(-3.14,3.14,3),'xticklabel',{'-\pi','0','\pi'})
        end
        
        gos = [gos V(rec).var.hit{var} V(rec).var.miss{var}];
        nogos = [nogos V(rec).var.CR{var} V(rec).var.FA{var}];
        
        gthetas{rec}=tmpgos;
        ngthetas{rec}=tmpnogos;
        

    end
    gothetas = histc([gos],range);
    nogothetas = histc([nogos],range);
    figure(89+d); clf; 
    bar(range,gothetas,'b');alpha(.5)
    hold on;bar(range,nogothetas,'r');alpha(.5)
    if var == 5
        hold on; bar(range,nogothetas,'r');alpha(.5)
        set(gca,'xlim',[min(range) max(range)],'xtick',linspace(-3.14,3.14,5),'xticklabel',{'-\pi','-\pi/2','0','\pi/2','\pi'})
    end
    ylabel('Number of Touches')
%     print(figure(31),'-dtiff',['Z:\Users\Jon\Projects\Characterization\' U{rec}.meta.layer '\Figures\'  U{rec}.meta.layer '_' V(1).varNames{var} 'DistPOP'])
%     print(figure(32),'-dtiff',['Z:\Users\Jon\Projects\Characterization\' U{rec}.meta.layer '\Figures\'  U{rec}.meta.layer '_' V(1).varNames{var} 'DistINDIV'])
% %     close all


end


thetameandiff = cellfun(@mean,ngthetas)-cellfun(@mean,gthetas);


gowidth = groupgorange(:,2)-groupgorange(:,1);
nogowidth = groupnogorange(:,2)-groupnogorange(:,1);


figure(321);hold on;h=scatter(iv{d},thetameandiff,'filled');
h.CData = rgb(colors{d});

thetacollat = [thetacollat thetameandiff];

end

%Plotting for SBIASxtouch count distribution figure320
% Plotting for SBIAS x meanof theta distributoin figure 321
xdata = cell2mat(iv');
ydata = collat';
thetaydata = thetacollat';
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
legend('Discrete','Semi-Continuous','Continuous','location','northwest')
 set(gca,'ylim',[-1 1],'ytick',[-1:.5:1],'xlim',[-5 25])
disp(['Rsquared = ' num2str(adjrsq) ' and pvalue = ' num2str(pval)])
xlabel('Go NoGo Gap (theta)')
ylabel('More nogo touches ------ More go touches')
set(figure(320), 'Units', 'pixels', 'Position', [0, 0, 600, 750]);

figure(321); plot(xdata(orders),t,'k');
legend('Discrete','Semi-Continuous','Continuous','location','northwest')
 set(gca,'xlim',[-5 30])
disp(['Rsquared = ' num2str(tadjrsq) ' and pvalue = ' num2str(tpval)])
xlabel('Gap (theta)')
ylabel('Difference between mean of go - mean of nogo')
set(figure(321), 'Units', 'pixels', 'Position', [0, 0, 600, 750]);


%% Task D x SBIAS CORRELATIOn 
ydata =  cell2mat(POP.SBIAS')
xdata = cell2mat(POP.taskD')


[~, orders] = sort(xdata);

[coeff, ~ , mu] = polyfit(xdata(orders),ydata(orders),1);
g = polyval(coeff,xdata(orders),[],mu);

modelvals = fitlm(xdata,ydata);
pval = modelvals.Coefficients{2,4};
adjrsq = modelvals.Rsquared.Adjusted;
disp(['Rsquared = ' num2str(adjrsq) ' and pvalue = ' num2str(pval)])

colors = {'DarkGreen','DarkMagenta','DarkTurquoise'};
ranges = [1:10;11:20;21:30]
figure(23);clf;
for d = 1:length(colors)
hold on; h=scatter(xdata(ranges(d,:)),ydata(ranges(d,:)),'filled');
h.CData = rgb(colors{d});
end
hold on; plot(xdata(orders),g,'k')
xlabel('Gap between go and nogo (theta)');ylabel('Search Bias')

legend('Discrete','Semi-Continuous','Continuous','location','northwest')


