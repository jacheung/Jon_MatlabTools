

%% AIC
%use exp(a0)*Delta to calculate FR/ms as the results are in log. exp gets rid of

type = {SM,BV};
for d = 2
    U=type{d} ;
    clear V
    [V] = classifierWrapper(U);
    
%     uberedwt = POPv2.logfeatures.ubered{d};
%     countswt = POPv2.logfeatures.counts{d};
%     thetawt = POPv2.logfeatures.theta{d};
    uberedwt = uberedfeat{d};
    thetawt = thetafeat{d};
    countswt = countsfeat{d};
    timingwt = timingfeat{d};
    motorwt = motorfeat{d};

    classes = 'lick';
    % 1) 'gonogo' 2) 'FAvsCR' 3) 'lick' 4) allBehavTypes
    
    sample ='bias';
    % 1) 'bias' (takes 70% from each class for train) 2) 'random' just takes
    % random 70% to train
    
    % Only for 'ubered' or 'pas'
    normalization = 'whiten';
    % 1) 'whiten' 2) 'none';
    
    % Only for 'ubered'
    removal = 'no';
    % 1) 'yes' to remove 0 touch trials and auto classify as CR
    
    balance = 'off';
    
    nanMethod = 'random';
    % 1) random (resample NaN vars from all touches)
    % 2) peakPro (replace NaN using var from max protraction)
    % 3) resampSameD (replace NaN using vas from touches in same trial type)
    dropempty = 'yes';
    % 1) 'yes' = drop trials with 0 touches
    % 2) 'no' = keep all trials
    
    for k = 1:size(uberedwt,2)
        designvars = 'ubered';
        % 1) 'theta' 2) 'pas' (phase amp midpoint) 3) 'counts' $) 'ubered'
        mouseNum = k;
        
        [DmatX, DmatY, motorX] = designMatrixBuilderv2(V(mouseNum ),U{mouseNum },designvars,classes,normalization,removal,balance,nanMethod,dropempty);
        thetavar = DmatX(:,1);
        countsvar = DmatX(:,2);
        yvals = double((DmatY==1));
        
        %         uberedloglikelihood =sum( -log(1+exp(uberedwt(3,mouseNum) + uberedwt(2,mouseNum)*thetavar + uberedwt(3,mouseNum)*countsvar)) + (yvals .* (uberedwt(3,mouseNum) + uberedwt(2,mouseNum)*thetavar + uberedwt(3,mouseNum)*countsvar)));
        uberedloglikelihood =sum( -log(1+exp(uberedwt(3,mouseNum) + uberedwt(2,mouseNum)*thetavar + uberedwt(3,mouseNum)*countsvar))) + sum(yvals .* (uberedwt(3,mouseNum) + uberedwt(2,mouseNum)*thetavar + uberedwt(3,mouseNum)*countsvar));
        uberedparams = size(DmatX,2);
        uberedsamples = size(DmatX,1);
        
        uberedaic{d}(k)=(-2*uberedloglikelihood)+(2*uberedparams);
        uberedbic{d}(k)=(-2*uberedloglikelihood)+(uberedparams*log(uberedsamples));
    end
    
    for k = 1:size(countswt,2)
        designvars = 'counts';
        % 1) 'theta' 2) 'pas' (phase amp midpoint) 3) 'counts' $) 'ubered'
        mouseNum = k;
        
        [DmatX, DmatY, motorX] = designMatrixBuilderv2(V(mouseNum ),U{mouseNum },designvars,classes,normalization,removal,balance,nanMethod,dropempty);
        xvar = DmatX(:,1);
        yvals = double((DmatY==1));
        
        %       countsloglikelihood =sum( -log(1+ exp(countswt(2,mouseNum) + countswt(1,mouseNum)*xvar) + (yvals .* (countswt(2,mouseNum) + countswt(1,mouseNum)*xvar)) ));
        countsloglikelihood =sum( -log(1+ exp(countswt(2,mouseNum) + countswt(1,mouseNum)*xvar))) + sum(yvals .* (countswt(2,mouseNum) + countswt(1,mouseNum)*xvar));
        countsparams = size(DmatX,2);
        
        xsamples = size(DmatX,1);
        countsaic{d}(k)=(-2* countsloglikelihood)+(2* countsparams);
        countsbic{d}(k)=(-2* countsloglikelihood)+(countsparams*log(xsamples));
    end
    
    for k = 1:size(thetawt,2)
        designvars = 'theta';
        % 1) 'theta' 2) 'pas' (phase amp midpoint) 3) 'counts' $) 'ubered'
        mouseNum = k;
        
        [DmatX, DmatY, motorX] = designMatrixBuilderv2(V(mouseNum ),U{mouseNum },designvars,classes,normalization,removal,balance,nanMethod,dropempty);
        xvar = DmatX(:,1);
        yvals = double((DmatY==1));
        
        %       thetaloglikelihood =sum( -log(1+exp(thetawt(2,mouseNum) + thetawt(1,mouseNum)*xvar) + (yvals .* (thetawt(2,mouseNum) + thetawt(1,mouseNum)*xvar))));
        thetaloglikelihood =sum( -log(1+ exp(thetawt(2,mouseNum) + thetawt(1,mouseNum)*xvar))) + sum(yvals .* (thetawt(2,mouseNum) + thetawt(1,mouseNum)*xvar));
        thetaparams = size(DmatX,2);
        
        
        xsamples = size(DmatX,1);
        thetaaic{d}(k)=(-2* thetaloglikelihood)+(2* thetaparams);
        thetabic{d}(k)=(-2* thetaloglikelihood)+(thetaparams*log(xsamples));
    end
    
    for k = 1:size(timingwt,2)
        designvars = 'timing';
        % 1) 'theta' 2) 'pas' (phase amp midpoint) 3) 'counts' $) 'ubered'
        mouseNum = k;
        
        [DmatX, DmatY, motorX] = designMatrixBuilderv2(V(mouseNum ),U{mouseNum },designvars,classes,normalization,removal,balance,nanMethod,dropempty);
        xvar = DmatX(:,1);
        yvals = double((DmatY==1));
        
        %       thetaloglikelihood =sum( -log(1+exp(thetawt(2,mouseNum) + thetawt(1,mouseNum)*xvar) + (yvals .* (thetawt(2,mouseNum) + thetawt(1,mouseNum)*xvar))));
        thetaloglikelihood =sum( -log(1+ exp(timingwt(2,mouseNum) + timingwt(1,mouseNum)*xvar))) + sum(yvals .* (timingwt(2,mouseNum) + timingwt(1,mouseNum)*xvar));
        thetaparams = size(DmatX,2);
        
        
        xsamples = size(DmatX,1);
        timingaic{d}(k)=(-2* thetaloglikelihood)+(2* thetaparams);
        timingbic{d}(k)=(-2* thetaloglikelihood)+(thetaparams*log(xsamples));
    end
    
    for k = 1:size(motorwt,2)
        designvars = 'motor';
        % 1) 'theta' 2) 'pas' (phase amp midpoint) 3) 'counts' $) 'ubered'
        mouseNum = k;
        
        [DmatX, DmatY, motorX] = designMatrixBuilderv2(V(mouseNum ),U{mouseNum },designvars,classes,normalization,removal,balance,nanMethod,dropempty);
        xvar = DmatX(:,1);
        yvals = double((DmatY==1));
        
        %       thetaloglikelihood =sum( -log(1+exp(thetawt(2,mouseNum) + thetawt(1,mouseNum)*xvar) + (yvals .* (thetawt(2,mouseNum) + thetawt(1,mouseNum)*xvar))));
        thetaloglikelihood =sum( -log(1+ exp(motorwt(2,mouseNum) + motorwt(1,mouseNum)*xvar))) + sum(yvals .* (motorwt(2,mouseNum) + motorwt(1,mouseNum)*xvar));
        thetaparams = size(DmatX,2);
        
        
        xsamples = size(DmatX,1);
        motoraic{d}(k)=(-2* thetaloglikelihood)+(2* thetaparams);
        motorbic{d}(k)=(-2* thetaloglikelihood)+(thetaparams*log(xsamples));
    end
    
    
end

% semiaic = [countsaic{1}; thetaaic{1}; timingaic{1}; motoraic{1} ;uberedaic{1}];
contaic = [countsaic{2}; thetaaic{2}; timingaic{2};motoraic{2} ;uberedaic{2}];
% semibic = [countsbic{1}; thetabic{1}; timingbic{1};uberedbic{1}];
% contbic = [countsbic{2}; thetabic{2}; timingbic{2};uberedbic{2}];

% normsemi = (semiaic(1:5,:)./max(semiaic(1:5,:)))';
normcont = (contaic(1:5,:)./max(contaic(1:5,:)))';

% znormsemi = (semiaic(1:4,:)./max(semiaic(1:4,:)))';
znormcont = (contaic(1:4,:)./max(contaic(1:4,:)))';
% znormsemi(znormsemi==0)=nan;
[~,sidx] = min(semiaic);
[~,cidx] = min(contaic);

colors = {'Gold','DarkTurquoise'};
figure(465);clf;
subplot(1,2,1)
hold on;
plot(1:5,normsemi','o','Markerfacecolor',rgb(colors{1}),'Markeredgecolor',rgb(colors{1}),'color',rgb(colors{1}));
set(gca,'xtick',[],'ytick',0:.5:1,'ylim',[0 1])
% plot(1:5,normcont,'o','Markerfacecolor',rgb(colors{2}),'Markeredgecolor',rgb(colors{2}),'color',rgb(colors{2}));
set(gca,'xlim',[.5 5.5],'xtick',[],'ytick',0:.5:1,'ylim',[0 1])
axis square
% errorbar(1:5,nanmean(normcont),nanstd(normcont),'ko','markerfacecolor',rgb(colors{2}))
errorbar(1:5,nanmean(normsemi),nanstd(normsemi),'ko','markerfacecolor',rgb(colors{1}))



subplot(1,2,2)
hold on;
plot(1:4,znormsemi','o','Markerfacecolor',rgb(colors{1}),'Markeredgecolor',rgb(colors{1}),'color',rgb(colors{1}));
set(gca,'xtick',[],'ytick',0:.5:1,'ylim',[0 1])
% plot(1:4,znormcont,'o','Markerfacecolor',rgb(colors{2}),'Markeredgecolor',rgb(colors{2}),'color',rgb(colors{2}));
set(gca,'xlim',[.5 4.5],'xtick',[],'ytick',0:.5:1,'ylim',[0 1])
axis square
% errorbar(1:4,nanmean(znormcont),nanstd(znormcont),'ko','markerfacecolor',rgb(colors{2}))
errorbar(1:4,nanmean(znormsemi),nanstd(znormsemi),'ko','markerfacecolor',rgb(colors{1}))

figure(10);clf
b=bar([histcounts(sidx,1:5)./numel(sidx);histcounts(cidx,1:5)./numel(cidx)]);
set(gca,'ylim',[0 1],'ytick',0:.5:1,'xtick',[])

figure(11);clf
errorbar(1:4,nanmean(znormcont),nanstd(znormcont),'ko')
set(gca,'xlim',[0 5])


