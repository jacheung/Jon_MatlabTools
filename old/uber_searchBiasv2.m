ranges = -75:1:75;
% sbiastwo = nan(length(U),2);
% sbiasthree = nan(length(U),1);
% sbiasfour = nan(length(U),1);

type = {SM,BV};
% type = {Nx,BVx};
clearvars -except U D BV SM BVx Nx type ranges POP
close all
for p=1:length(type)   
U = type{p};


%%   
    for rec = 1:length(U)
        
        %choosing trials with which to evaluate exploration
        if strcmp(U{rec}.meta.layer,'D')
            farthestnogoT = find(U{rec}.meta.trialType==0);
        else
            farthestnogoT= find(U{rec}.meta.motorPosition<min(U{rec}.meta.motorPosition)+10000);
        end
        thetasall = squeeze(U{rec}.S_ctk(1,:,:));
        thetasnan = nan(size(thetasall));
        thetasnan(:,farthestnogoT) = 1;
        thetas = thetasall.*thetasnan;
        
        [whisks] = findMaxMinProtraction(U{rec},'sampling');
        
        selectedT=thetas(whisks.peakidx);
        selectedT(isnan(selectedT))=[];
        selectedT=sort(selectedT);
        [v,selectedT] = ecdf(selectedT);
        
        %moments of exploration distribution 
        means{p}(rec) = mean(selectedT);
        stds{p}(rec) = std(selectedT);
        skewness{p}(rec)= 3*(mean(selectedT)-median(selectedT))/std(selectedT);
        
        
        %Plots for looking at consistency of peak theta values across all
        %selected trials
%         thetasunrolled = thetasall(:);
%         nansunrolled = nan(length(thetasall(:)),1);
%         nansunrolled(U{rec}.peakIdx)=1;
%         finalthetas = thetasunrolled.*nansunrolled;
%         reshapedthetas = reshape(finalthetas,U{rec}.t,U{rec}.k)';
%         figure(5050);clf;imagesc(reshapedthetas(farthestnogoT,:));
%         colorbar
%         figure(4040);clf;errorbar(nanmean(reshapedthetas(farthestnogoT,:),2),1:size(reshapedthetas(farthestnogoT,:),1),nanstd(reshapedthetas(farthestnogoT,:),[],2),'ko','horizontal')
%         
        %peak whisk distribution
        [histovals ] = binslin(selectedT,selectedT,'equalE',numel(ranges),ranges(1),ranges(end));
        peakwhisk = cellfun(@numel,histovals)./ sum(cellfun(@numel,histovals));
        %p(whisk>angle) distribution
        [sorted, ~ , ~]=binslin(selectedT,v,'equalE',numel(ranges),ranges(1),ranges(end));
        pwhisk = 1-cellfun(@mean,sorted);
        
        [~,maxidx] = max(pwhisk);
        [~,minidx] = min(pwhisk);
        pwhisk(maxidx:minidx) = smooth(pwhisk(maxidx:minidx),5);
        pwhisk(pwhisk>1) = 1; %constraining since smooth sometimes makes vals above 1
        
        %p(pole at theta) distribution
        ppoles = nan(U{rec}.k,2);
        for b = 1:U{rec}.k
            touchIdx = [find(U{rec}.S_ctk(9,:,b)==1) find(U{rec}.S_ctk(12,:,b)==1)];
            if ~isempty(touchIdx)
                tmpThetas = U{rec}.S_ctk(1,touchIdx,b);
                ppoles(b,:) = [median(tmpThetas) std(tmpThetas)];
            end
        end
        
        %finding proportion of nogo trials that have at least one contact
        pnogotouches{p}(rec) = 1 - numel(intersect(find(isnan(ppoles(:,1))==1),find(U{rec}.meta.trialType == 0))) / numel(find(U{rec}.meta.trialType == 0));
        
        %some nogo trials don't have touches but need potential theta. Doing 2nd
        %order polynomial fit and then grabbing that value for expected theta
            % FOR DISCRETE
        if strcmp(U{rec}.meta.layer,'D')
            goavgtheta = mean(ppoles(U{rec}.meta.trialType==1,1));
            nogoavgtheta = nanmean(ppoles(U{rec}.meta.trialType==0,1));
            missedGoIdx = logical(isnan(ppoles(:,1))'.*double(U{rec}.meta.trialType == 1));
            missedNogoIdx = logical(isnan(ppoles(:,1))'.*double(U{rec}.meta.trialType == 0));
            ppoles(missedNogoIdx,1) = repmat(nogoavgtheta,sum(missedNogoIdx),1);
            ppoles(missedGoIdx,1) = repmat(goavgtheta,sum(missedGoIdx),1);
            
        elseif strcmp(U{rec}.meta.layer,'BV') || strcmp(U{rec}.meta.layer,'SM') || strcmp(U{rec}.meta.layer,'BVx')
            
            % FOR CONTINUOUS/SEMI
            polyinputs = sortrows([U{rec}.meta.motorPosition'  ppoles(:,1)]);
            polyinputs(isnan(polyinputs(:,2)),:)=[];
            [coeff, ~ , mu] = polyfit(polyinputs(:,1),polyinputs(:,2),2);
            thetafills = polyval(coeff,sort(U{rec}.meta.motorPosition),[],mu);
            
            % checking pole pos and fitted line
            % figure(650);clf;scatter(U{rec}.meta.motorPosition,ppoles(:,1))
            % hold on; plot(sort(U{rec}.meta.motorPosition),thetafills,'b')
            missedMotorThetas = U{rec}.meta.motorPosition(isnan(ppoles(:,1)));%find idx of nan vals
            filledMotorThetas = polyval(coeff,missedMotorThetas,[],mu); %find theta for motor idx
            ppoles(isnan(ppoles(:,1)))=filledMotorThetas; %fill nan vals
        end
        %bin pole distribution based on theta by 1 degree and transform to
        %probability
        [sorted2, sortedBy ,binBounds]=binslin(ppoles(:,1),ppoles(:,1),'equalE',numel(ranges),ranges(1),ranges(end));
        ppoletheta = cellfun(@numel,sorted2)./sum(cellfun(@numel,sorted2));
        
        gopoles = ppoles(U{rec}.meta.trialType==1,1);
        nogopoles = ppoles(U{rec}.meta.trialType==0,1);
        gos = round([nanmedian(gopoles)-1.5*nanstd(gopoles)  nanmedian(gopoles)+1.5*nanstd(gopoles) ])
        nogos = round([nanmedian(nogopoles)-1.5*nanstd(nogopoles)  nanmedian(nogopoles)+1.5*nanstd(nogopoles) ])
        
%         gos = round([min(gopoles) max(gopoles)])
%         nogos = round([min(nogopoles) max(nogopoles)])

        gos = positions{rec}(1,:);
        nogos = positions{rec}(2,:);
        
        goidx = gos(1):1:gos(end);
        nogoidx = nogos(1):1:nogos(end);
        
        
        %calculating "joint probability"
        jointProb = pwhisk.*ppoletheta;
        jointProb(isnan(jointProb))=0;
        
        [~,di] = intersect(ranges,goidx);
        gostuff{p}{rec} = [goidx' pwhisk(di) ppoletheta(di) jointProb(di)];        
        [~,do] = intersect(ranges,nogoidx);
        nogostuff{p}{rec} = [nogoidx' pwhisk(do) ppoletheta(do) jointProb(do)];
        

        
        %finding proportion of whisks that go past dbound
        dbound = floor(mean([di(end) do(1)]));
        whisk{p}(rec) = pwhisk(dbound);
        %proportion of whisks that go past closest nogo
%         whisk{p}(rec) = nogostuff{p}{rec}(1,2);
        
        %sbias raw variables
        sbiastwo{p}(rec) = (nansum(gostuff{p}{rec}(:,4)) - nansum(nogostuff{p}{rec}(:,4)) ) ./(nansum(gostuff{p}{rec}(:,4)) + nansum(nogostuff{p}{rec}(:,4)) );
        sbiasequal{p}(rec) = 1 - nanmean(nogostuff{p}{rec}(:,2));
        sbiasthree{p}(rec) = nansum(nogostuff{p}{rec}(:,4)) ./ (nansum(nogostuff{p}{rec}(:,4)) + nansum(gostuff{p}{rec}(:,4))); %proportion of nogo/total contact
        sbiasfour{p}(rec) = 1-(nansum(nogostuff{p}{rec}(:,4))./nansum(nogostuff{p}{rec}(:,3))); %proportion of nogo contacted
    end
end
%% plotting sbias across different tasks
colors = {'DarkGreen','DarkMagenta','DarkTurquoise'};
figure(4050);clf
figure(54020);clf
feature = sbiasequal;
for d = 1:3
    figure(4050);hold on;
    scatter(feature{d},ones(length(feature{d}),1).*d,'filled','markerfacecolor',rgb(colors{d}));
    hold on; errorbar(mean(feature{d}),d,std(feature{d}),'horizontal','ko','markersize',15,'markerfacecolor',rgb(colors{d}));

    %checking which variable if whisk past db is correlated with proportion
    %of nogo trials being contacted 
%     figure(54020); hold on
%     scatter(gngtaskease{d},gngcountsratios{d},'filled','markerfacecolor',rgb(colors{d}));
    
end

get(figure(4050))
set( gca,'ydir','reverse','ylim',[.5 3.5],'ytick',[],'xlim',[0 1],'xtick',[0:.5:1])
set(figure(4050), 'Units', 'pixels', 'Position', [0, 0, 500, 650]);
xlabel('Proportion of nogo trials with a contact');


get(figure(54020));
set(gca,'xtick',[-10:10:40],'xlim',[-10 40],'ylim',[-1 1],'ytick',[-1 0 1])
% plot([0 1],[0 0],'-.k')
xlabel('Task ease')
ylabel('count ratios');
set(figure(54020), 'Units', 'pixels', 'Position', [0, 0, 500, 650]);

xax = cell2mat(gngtaskease)';
ydata = cell2mat(gngcountsratios)';
[~, orders] = sort(xax);

[coeff, ~ , mu] = polyfit(xax(orders),ydata(orders),1);
f = polyval(coeff,xax(orders),[],mu);
hold on; plot(xax(orders),f,'k');

modelvals = fitlm(xax,ydata);
modelvals.Coefficients(2,4) ;
pval = modelvals.Coefficients{2,4}
adjrsq = modelvals.Rsquared.Adjusted


%%
colors = {'DarkGreen','DarkMagenta','DarkTurquoise'};
popfeatdom = [];
figure(38);clf
%SBIAS testing
feature = sbiasfour;
for d = 1:3
    
    
%         xax = POP.SBIAS{d}; %OG
    % xax = sbiasthree;
    xax = feature{d};
%     xax = pnogotouches{d};
    
    THETAS = POP.feature{d}(:,1);
    COUNTS = POP.feature{d}(:,2);
    featureDom = (abs(THETAS)-abs(COUNTS))./(abs(THETAS)+abs(COUNTS));
    
    figure(38);
    hold on; h = scatter(xax,featureDom,'filled');
    h.CData = rgb(colors{d});
    
    popfeatdom = [popfeatdom ;featureDom];
    
end

xax = cell2mat(feature)';
% xax = cell2mat(POP.SBIAS');
% ydata = featureDom;
ydata = popfeatdom
[~, orders] = sort(xax);

[coeff, ~ , mu] = polyfit(xax(orders),ydata(orders),1);
f = polyval(coeff,xax(orders),[],mu);
hold on; plot(xax(orders),f,'k');
hold on; plot([-20 20],[0 0],'-.k')
set(gca,'ylim',[-1 1],'xlim',[0 1],'xtick',[-.5:.5:.5],'ytick',[-1:.5:1])
set(figure(38), 'Units', 'pixels', 'Position', [0, 0, 500, 650]);

modelvals = fitlm(xax,ydata);
modelvals.Coefficients(2,4) ;
% legend('Discrete','Semi-Continuous','Continuous')
pval = modelvals.Coefficients{2,4}
adjrsq = modelvals.Rsquared.Adjusted
%%








%Plotting
scaleval = .05;
normalizedpwhisk = pwhisk.*scaleval;

figure(5080);clf
bar(1:length(peakwhisk),peakwhisk,'facecolor',[.8 .8 .8]);
set(gca,'ylim',[0 scaleval],'xtick',1:25:150,'xticklabel',ranges(1):25:ranges(end),'xlim',[50 125])
ylabel('Proportion of whisks');xlabel('Theta')
hold on; plot(1:length(pwhisk),normalizedpwhisk,'k');
a2 = axes('YAxisLocation', 'right');
set(a2,'color','none');set(a2,'XTick',[]);
set(gca,'YLim',[0 1],'ytick',[0:.25:1]);ylabel('Probability')

figure(504380);clf
%  %plotting eCDF of whisks>angle
[sorted3, sortedby ,~]=binslin(di,ppoletheta(di),'equalE',6,di(1),di(end));%plotting p(pole go)
hold on; bar(cellfun(@mean,sortedby),cellfun(@mean,sorted3),'b')


[sorted4, sortedby4 ,~]=binslin(do,ppoletheta(do),'equalE',6,do(1),do(end));%plotting p(pole nogo)
hold on; bar(cellfun(@mean,sortedby4),cellfun(@mean,sorted4),'r')


[sorted4, sortedby4 ,~]=binslin(do,ppoletheta(do),'equalE',6,do(1),do(end));%plotting p(pole nogo)

% hold on; bar(di(1):di(end),ppoletheta(di),'b'); %plotting p(pole go)
% hold on; bar(do(1):do(end),ppoletheta(do),'r'); %plotting p(pole nogo)
hold on; bar(1:length(jointProb),jointProb,'g'); %plotting joint distribution P(whisk>theta,pole==theta)
set(gca,'xlim',[50 125],'ylim',[0 .1],'xtick',0:25:150,'xticklabel',ranges(1):25:ranges(end))
ylabel('Probability');xlabel('Theta')

plot(1:length(pwhisk),pwhisk.*.1,'k');
legend('P(go pole)','P(nogo pole)','Contact likelihood?','P(whisk>theta)','location','best')


a2 = axes('YAxisLocation', 'right');
set(a2,'color','none');set(a2,'XTick',[]);
set(gca,'YLim',[0 1],'ytick',[0:.25:1]);ylabel('Probability')







