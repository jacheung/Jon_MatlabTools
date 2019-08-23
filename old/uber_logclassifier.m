


clear V
clear F1G RMSEgroup
[V] = classifierWrapper(U);

%% PARAMETERS SETTING

clearvars -except V U BV D SM F1G RMSEgroup RMSEdecomp
numIterations = 10;

designvars = 'ubered';
% 1) 'theta' 2) 'pas' (phase amp midpoint) 3) 'counts' $) 'ubered'
classes = 'lick';
% 1) 'gonogo' 2) 'FAvsCR' 3) 'lick' 4) allBehavTypes

sample ='bias';
% 1) 'bias' (takes 70% from each class for train) 2) 'random' just takes
% random 70% to train

% Only for 'ubered' or 'pas'
normalization = 'none';
% 1) 'whiten' 2) 'none';

% Only for 'ubered'
removal = 'yes';
% 1) 'yes' to remove 0 touch trials and auto classify as CR

balance = 'off';
% Design Matrix


accprop=cell(1,length(V));

for rec = 1:length(U)
    
    if strcmp(U{rec}.meta.layer,'D')
        U{rec}.meta.motorPosition = U{rec}.meta.trialType;
        U{rec}.meta.motorPosition(U{rec}.meta.motorPosition == 0) = -1;
    end
    motorPos = U{rec}.meta.motorPosition;
    hxMotor = motorPos(logical(V(rec).trialNums.matrix(1,:)))';
    FAxMotor = motorPos(logical(V(rec).trialNums.matrix(3,:)))';
    CRxMotor = motorPos(logical(V(rec).trialNums.matrix(4,:)))';
    %
    switch designvars
        case 'theta'
            hx = [V(rec).var.hit{1}'];
            FAx = [V(rec).var.FA{1}'];
            CRx = [V(rec).var.CR{1}'];
            hy = ones(size(hx,1),1);
            FAy = ones(size(FAx,1),1)+1;
            CRy1 = ones(size(CRx,1),1)+1;
            
            
            htmp =[V(rec).touchNum.hit' hxMotor];
            FAtmp = [V(rec).touchNum.FA' FAxMotor];
            CRtmp = [V(rec).touchNum.CR' CRxMotor];
            
            %             htmp =[V(rec).touchNum.hit' ];
            %             FAtmp = [V(rec).touchNum.FA' ];
            %             CRtmp = [V(rec).touchNum.CR' ];
            %
            newMvals = {htmp,FAtmp,CRtmp};
            for d = 1:length(newMvals)%for each trial type
                newmotor = [] ;
                for k=1:length(newMvals{d}) %for trials...
                    j = newMvals{d}(k); %for the number of touches within trial
                    if ~j==0
                        tmp = repmat(newMvals{d}(k,2),j,1);
                        newmotor = [newmotor; tmp];
                    end
                end
                if d == 1
                    hx = [hx newmotor];
                elseif d == 2
                    FAx = [FAx newmotor];
                elseif d == 3
                    CRx = [CRx newmotor];
                end
            end
            
            
        case 'pas'
            ntmp=find(V(rec).var.hit{5}<=0);ptmp=find(V(rec).var.hit{5}>0);
            hx = [V(rec).var.hit{3}(ntmp)' V(rec).var.hit{4}(ntmp)' V(rec).var.hit{5}(ntmp)'];
            hy = ones(size(hx,1),1);
            Fntmp=find(V(rec).var.FA{5}<=0);Fptmp=find(V(rec).var.FA{5}>0);
            FAx = [V(rec).var.FA{3}(Fntmp)' V(rec).var.FA{4}(Fntmp)' V(rec).var.FA{5}(Fntmp)'];
            FAy = ones(size(FAx,1),1)+1;
            Cntmp=find(V(rec).var.CR{5}<=0);Cptmp=find(V(rec).var.CR{5}>0);
            CRx = [V(rec).var.CR{3}(Cntmp)' V(rec).var.CR{4}(Cntmp)' V(rec).var.CR{5}(Cntmp)'];
            CRy1 = ones(size(CRx,1),1)+1;
            %             hx = [hx hxMotor];
            %             FAx = [FAx FAxMotor];
            %             CRx = [CRx CRxMotor] ;
        case 'counts'
            
            hx = [V(rec).touchNum.hit'];
            hy = ones(size(hx,1),1);
            FAx = [V(rec).touchNum.FA'];
            FAy = ones(size(FAx,1),1)+1;
            CRx = [V(rec).touchNum.CR'];
            CRy1 = ones(size(CRx,1),1)+1;
            hx = [hx hxMotor];
            FAx = [FAx FAxMotor];
            CRx = [CRx CRxMotor] ;
            
            
        case 'ubered'
            
            [hx,mx,FAx,CRx] = meanVarfinder (V(rec),1,'NaN');
            hx = [hx' V(rec).touchNum.hit' V(rec).licks.oneT.hit' V(rec).licks.oneT.hit'.*V(rec).licks.twoT.hit' V(rec).licks.oneT.hit'.*V(rec).licks.twoT.hit'.*V(rec).licks.threeT.hit'];
            %                 mx = [mx' V(rec).touchNum.miss' V(rec).licks.oneT.miss'];
            FAx = [FAx' V(rec).touchNum.FA' V(rec).licks.oneT.FA' V(rec).licks.oneT.FA'.*V(rec).licks.twoT.FA' V(rec).licks.oneT.FA'.*V(rec).licks.twoT.FA'.*V(rec).licks.threeT.FA'];
            CRx = [CRx' V(rec).touchNum.CR' V(rec).licks.oneT.CR' V(rec).licks.oneT.CR'.*V(rec).licks.twoT.CR' V(rec).licks.oneT.CR'.*V(rec).licks.twoT.CR'.*V(rec).licks.threeT.CR'];
            
            hy = ones(size(hx,1),1);
            %                 my = ones(size(mx,1),1);
            FAy = ones(size(FAx,1),1)+1;
            CRy1 = ones(size(CRx,1),1)+1;
            
            hx = [hx hxMotor];
            FAx = [FAx FAxMotor];
            CRx = [CRx CRxMotor] ;
            
            switch removal
                case 'on'
                    find(DmatX(:,2) == 0)
                    
            end
    end
    
    
    switch classes
        case 'gonogo'
            DmatX=[hx;FAx;CRx]; DmatY=[hy;FAy;CRy1]; %complete design matrix for govsnogo
            colors = {'b','r'};
            g1 = [hy]; g2 = [FAy ; CRy1];
        case 'FAvsCR'
            if strcmp(balance,'on')
                [FAx,FAy,CRx,CRy1] = FACRBalance(FAx,CRx);
                DmatX = [FAx;CRx]; DmatY = [FAy-1;CRy1];
            else
                DmatX = [FAx;CRx]; DmatY = [FAy-1;CRy1];
            end
            colors = {'g','r'};
        case 'lick'
            if strcmp(balance,'on')
                [lix,lixy,nolixx,nolixy] = FACRBalance([hx ; FAx],CRx);
                DmatX = [lix(:,1:size(lix,2)-1);nolixx(:,1:size(lix,2)-1)]; DmatY = [lixy-1;nolixy];
                motorX = [lix(:,size(lix,2));nolixx(:,size(lix,2))];
            else
                DmatX = [hx(:,1:size(hx,2)-1);FAx(:,1:size(FAx,2)-1);CRx(:,1:size(CRx,2)-1)]; DmatY = [hy;FAy-1;CRy1];
                motorX = [hx(:,size(hx,2));FAx(:,size(FAx,2));CRx(:,size(CRx,2))];
            end
           
            
            if strcmp(designvars,'ubered') || strcmp(designvars,'pas')
                
                switch removal
                    case 'yes'
                        notouchidx= find(DmatX(:,2) == 0);
                        DmatX(notouchidx,:) = [];
                        DmatY(notouchidx,:) = [];
                end
                
                switch normalization
                    case 'whiten'
                        DmatX = filex_whiten(DmatX);
                end
                
                
            end
            
            g1= DmatY(DmatY == 1);
            g2= DmatY(DmatY == 2);
            
            colors = {'b','r'};
        case 'allBehavTypes'
            DmatX=[hx;FAx;CRx]; DmatY = [hy;FAy;CRy1+1];
    end
    
    [bagPred{rec} ,bagReal{rec},oobErrorBaggedEnsemble{rec}] = uber_baggingClassifier(DmatX,DmatY);
    
    clear Acc
    clear opt_thresh
    motorPlick = [];
    
    for reit = 1:numIterations
        rando = randperm(length(DmatX));
        tmpDmatX=DmatX(rando,:);tmpDmatY=DmatY(rando,:);
        if strcmp(classes,'lick')
            tmpMotorX = motorX(rando,:);
        end
        switch sample
            case 'bias'
                %             %FOR FA VS CR
                %             sample evenly from FA and CR for training set
                g1counts = round(numel(g1)*.7);
                g2counts = round(numel(g2)*.7);
                g1s = find(tmpDmatY==unique(g1));
                g2s = find(tmpDmatY==unique(g2));
                train=[g2s(1:g2counts);g1s(1:g1counts)];
                test = [1:length(tmpDmatY)]';
                test(train)=[];
                
                [thetas,cost,~] = ML_oneVsAll(tmpDmatX(train,:),tmpDmatY(train,:),numel(unique(DmatY)),0);
                Bfeat{rec}.theta{reit}=thetas;
                
                [pred,opt_thresh(reit),prob]=ML_predictOneVsAll(thetas,tmpDmatX(test,:)...
                    ,tmpDmatY(test,:),'Max');
                Acc(reit)= mean(double(pred == tmpDmatY(test))) * 100;
                F1s(reit,:) = F1score(pred,tmpDmatY(test),2);
                
                if strcmp(classes,'lick')
                    motorPlick= [motorPlick;tmpMotorX(test) prob];
                end
                accprop{rec}=[accprop{rec} ; pred tmpDmatY(test)];
                
            case 'random'
                % use this for GO vs NOGO since it doesnt bias sampling
                [thetas,cost,~] = ML_oneVsAll(tmpDmatX(1:end*.7,:),tmpDmatY(1:end*.7,:),numel(unique(DmatY)),0);
                Bfeat{rec}.theta{reit}=thetas;
                [pred,opt_thresh(reit),prob]=ML_predictOneVsAll(thetas,tmpDmatX(end*.7:end,:)...
                    ,tmpDmatY(end*.7:end,:),'Max');
                Acc(reit)= mean(double(pred == tmpDmatY(end*.7:end) )) * 100;
                F1s(reit,:) = F1score(pred,tmpDmatY(end*.7:end),2);
                
                if strcmp(classes,'lick')
                    motorPlick= [motorPlick;tmpMotorX(end*.7:end) prob];
                end
                accprop{rec}=[accprop{rec} ; pred tmpDmatY(end*.7:end)];
                
        end
        %fprintf('\nTraining Set Accuracy: %f\n', mean(double(pred == tmpDmatY(end*.6:end))) * 100)
        
    end
    if strcmp(classes,'lick')
        train_motorPlick{rec} = motorPlick;
    end
    train_F1s(rec,:) = nanmean(nansum(F1s,2));
    trainF1sstd(rec,:)=nanstd(nansum(F1s,2));
    train_predOpt(rec)=mean(opt_thresh);
    train_Acc(rec) = mean(Acc);
    train_std(rec) = std(Acc);


end


for rec = 1:length(V)
    figure(9);
    subplot(2,5,rec);plot(1:length(oobErrorBaggedEnsemble{rec}),oobErrorBaggedEnsemble{rec},'color',[.8 .8 .8])
    hold on; plot(1:length(oobErrorBaggedEnsemble{rec}),mean(oobErrorBaggedEnsemble{rec}),'k','linewidth',3)
    
    figure(10);hold on
    scatter(rec,mean(mean(bagPred{rec}==bagReal{rec})),'o')

end

if strcmp(designvars,'theta')
    F1G{1} = train_F1s ;
elseif strcmp(designvars,'counts')
    F1G{2} = train_F1s ;
    figure(20);hold on;f=scatter(F1G{1},F1G{2});
    if strcmp(U{rec}.meta.layer,'D')
        f.CData = [rgb('DarkGreen')];
    elseif strcmp(U{rec}.meta.layer,'SM')
        f.CData = [rgb('DarkMagenta')];
    elseif strcmp(U{rec}.meta.layer,'BV')
        f.CData = [rgb('DarkTurquoise')];
    end
    set(gca,'xlim',[0 2],'ylim',[0 2],'ytick',0:.5:2);
    xlabel('F1 Score Theta');ylabel('F1 Score Touch Count');
    %     hold on; plot([0 2],[0 2],'-.k')
    %     legend('Discrete','Semi Continuous','Continuous','location','northwest')
    
end


figure(22);clf;
hold on; errorbar(1:length(V),train_Acc,train_std,'ko')
xlabel('Mouse Number');ylabel('% Accuracy of Model')
set(gca,'ylim',[0 100])
xlim([.5 length(V)+.5])
%plotting TRUE POSITIVES (GUESSING class 1 and then it actually being class
%1)
predprop = zeros(length(V),3);
for recs = 1:length(V)
    guesses=accprop{recs};
    for i = 1:3
        tmp=find(guesses(:,2)==i);
        predprop(recs,i)=mean(guesses(tmp,1)==i)*100;
    end
end
hold on; scatter(1:length(V),predprop(:,1),'b');
scatter(1:length(V),predprop(:,2),'r')
legend('Model Accuracy','Lick Prediction Accuracy','No Lick Prediction Accuracy','location','southeast')
title([U{rec}.meta.layer ' ' designvars ' ' classes])
print(figure(22),'-dtiff',['Z:\Users\Jon\Projects\Characterization\' U{rec}.meta.layer '\Figures\'  U{rec}.meta.layer '_' classes '_' designvars '_0touch_removal_' removal])

%% FOR WEIGHTS
for rec = 1:length(V)
    ms=cell2mat(Bfeat{rec}.theta);
    optfeat(:,rec)=mean(reshape(ms(1,:)',size(DmatX,2)+1,numIterations),2);
    optfeatstd(:,rec)=std(reshape(ms(1,:)',size(DmatX,2)+1,numIterations),0,2);
end
figure(92);filex_barwitherr(optfeatstd',optfeat')

ylabel('Weight');
if strcmp(designvars,'pas')
    legend('Bias','Amplitude at Touch','Midpoint at Touch','Phase at Touch','Location','northeast')
elseif strcmp(designvars,'ubered')
    legend('Bias','Theta at Touch','Touch Count','Previous Trial Lick','2 Previous T Lick', '3 Previous T Lick','Location','northeast')
end

print(figure(92),'-dtiff',['Z:\Users\Jon\Projects\Characterization\' U{rec}.meta.layer '\Figures\'  U{rec}.meta.layer '_' classes '_' designvars '_WEIGHTS'])


%% Psychometric Curve Comparison b/t Model and Mouse
figure(3);clf
for rec = 1:length(V)
    
    if strcmp(U{rec}.meta.layer,'D')
        nolickmean = train_motorPlick{rec}(train_motorPlick{rec}(:,1) == -1,2);
        lickmean = train_motorPlick{rec}(train_motorPlick{rec}(:,1) == 1,2);
        gos = U{rec}.meta.trialType ==1;
        nogos = U{rec}.meta.trialType == 0 ;
        
        modelmean = [mean(nolickmean) mean(lickmean)];
        modelstd = [std(nolickmean) std(lickmean)];
        
        realmean = [1-mean(U{rec}.meta.trialCorrect(nogos)) mean(U{rec}.meta.trialCorrect(gos))];
        
        RMSE = sqrt(nanmean((realmean-modelmean).^2));
        
        figure(3);subplot(2,5,rec)
        scatter([-1 1],realmean,'filled','r')
        hold on;errorbar([-1 1],modelmean,modelstd,'ko')
        set(gca,'xlim',[-1.5 1.5],'ylim',[0 1])
        
        if rec ==5
            legend('Mouse Performance','Model Performance','location','southeast')
        elseif rec == 3
            title([U{rec}.meta.layer ' ' designvars ' ' classes])
        end
        
        text(.5,.4,['RMSE = ' num2str(RMSE)])
        
        
        
    elseif strcmp(U{1}.meta.layer,'SM')
        
        [sorted]= binslin(train_motorPlick{rec}(:,1),train_motorPlick{rec},'equalE',12,U{rec}.meta.ranges(1),U{rec}.meta.ranges(2));
        lickmean=cell2mat(cellfun(@mean,sorted,'uniformoutput',0));
        lickstd=cell2mat(cellfun(@std,sorted,'uniformoutput',0));
        xranges = U{rec}.meta.ranges(1):10000:U{rec}.meta.ranges(2);
        
        real=[U{rec}.meta.motorPosition;V(rec).trialNums.matrix(5,:)]';%only taking lick row and motorPos row
        [realsorted]= binslin(real(:,1),real,'equalE',12,U{rec}.meta.ranges(1),U{rec}.meta.ranges(2));
        reallickmean=cell2mat(cellfun(@(x) mean(x,1),realsorted,'uniformoutput',0));
        reallickstd=cell2mat(cellfun(@(x) std(x,0,1),realsorted,'uniformoutput',0));
        
        RMSE = sqrt(nanmean((reallickmean(:,2)-lickmean(:,2)).^2));
        
        
        figure(3);subplot(2,5,rec)
        xranges(isnan(lickmean(:,1)))=[];
        lickmean(isnan(lickmean(:,1)),:)=[];
        lickstd(isnan(lickstd(:,1)),:)=[];
        reallickmean(isnan(reallickmean(:,1)),:)=[];
        scatter(xranges(1),lickmean(1,2),'ko','filled')
        hold on;filex_shadedErrorBar(xranges(end-(length(lickmean)-2):end),lickmean(2:end,2),lickstd(2:end,2),'k')
        scatter(xranges(1),reallickmean(1,2),'ro','linewidth',1)
        plot(xranges(end-(length(lickmean)-2):end),reallickmean(2:end,2),'r','linewidth',2)
        
        if rec == 3
            title([U{rec}.meta.layer ' ' designvars ' ' classes])
        end
        
        set(gca,'xtick',[xranges(1) xranges(end-(length(lickmean)-2)) xranges(end)],'xticklabel',[-1 0 1],'ylim',[0 1],'xlim',[xranges(1)-5000 xranges(end)+5000])
        xlabel('Motor Position')
        ylabel('Lick Probability')
        text(xranges(end)*.75,.4,['RMSE = ' num2str(RMSE)])
        
    elseif strcmp(U{1}.meta.layer,'BV')
        [sorted]= binslin(train_motorPlick{rec}(:,1),train_motorPlick{rec},'equalE',12,U{rec}.meta.ranges(1),U{rec}.meta.ranges(2));
        lickmean=cell2mat(cellfun(@mean,sorted,'uniformoutput',0));
        lickstd=cell2mat(cellfun(@std,sorted,'uniformoutput',0));
        xranges = U{rec}.meta.ranges(1):10000:U{rec}.meta.ranges(2);
        
        real=[U{rec}.meta.motorPosition;V(rec).trialNums.matrix(5,:)]';%only taking lick row and motorPos row
        [realsorted]= binslin(real(:,1),real,'equalE',12,U{rec}.meta.ranges(1),U{rec}.meta.ranges(2));
        reallickmean=cell2mat(cellfun(@(x) mean(x,1),realsorted,'uniformoutput',0));
        reallickstd=cell2mat(cellfun(@(x) std(x,0,1),realsorted,'uniformoutput',0));
        RMSE = sqrt(nanmean((reallickmean(:,2)-lickmean(:,2)).^2));
        RMSEdec = sqrt((reallickmean(:,2)-lickmean(:,2)).^2);
        
        
        figure(3);subplot(2,5,rec)
        plot(xranges',reallickmean(:,2),'r','linewidth',1)
        hold on;filex_shadedErrorBar(xranges,lickmean(:,2),lickstd(:,2),'k');
        hold on; plot(xranges',reallickmean(:,2),'r','linewidth',5)
        set(gca,'xtick',[xranges(1) xranges(6) xranges(11)],'xticklabel',[-1 0 1],'ylim',[0 1],'xlim',[xranges(1)-5000 xranges(end)+5000])
        if rec ==5
            legend('Mouse Performance','Model Performance','location','southeast')
        elseif rec == 3
            title([U{rec}.meta.layer ' ' designvars ' ' classes])
            
        end
        
        xlabel('Motor Position')
        ylabel('Lick Probability')
        text(xranges(end)*.75,.4,['RMSE = ' num2str(RMSE)])
        
    end
    
    
    if strcmp(designvars,'counts')
        RMSEgroup(rec,1) = RMSE;
        RMSEdecomp{1}(:,rec) = RMSEdec;
    elseif strcmp(designvars,'theta')
        RMSEgroup(rec,2) = RMSE;
        RMSEdecomp{2}(:,rec) = RMSEdec;
    elseif strcmp(designvars,'ubered')
        RMSEgroup(rec,3) = RMSE;
        RMSEdecomp{3}(:,rec) = RMSEdec;
    end
    
    
    
end
set(gcf, 'Units', 'pixels', 'Position', [0, 0, 2000, 1000]);
print(figure(3),'-dtiff',['Z:\Users\Jon\Projects\Characterization\' U{rec}.meta.layer '\Figures\'  U{rec}.meta.layer '_' classes '_' designvars '_PSYCHO'])


%% Visualization of the decision boundaries.
params =1;
FCratio = [];
for rec = 1:length(U)
    ntmp=find(V(rec).var.hit{5}<=0);ptmp=find(V(rec).var.hit{5}>0);
    hx = [V(rec).var.hit{3}(ntmp)' V(rec).var.hit{4}(ntmp)' V(rec).var.hit{5}(ntmp)'];
    hy = ones(size(hx,1),1);
    Fntmp=find(V(rec).var.FA{5}<=0);Fptmp=find(V(rec).var.FA{5}>0);
    FAx = [V(rec).var.FA{3}(Fntmp)' V(rec).var.FA{4}(Fntmp)' V(rec).var.FA{5}(Fntmp)'];
    FAy = ones(size(FAx,1),1)+1;
    Cntmp=find(V(rec).var.CR{5}<=0);Cptmp=find(V(rec).var.CR{5}>0);
    CRx = [V(rec).var.CR{3}(Cntmp)' V(rec).var.CR{4}(Cntmp)' V(rec).var.CR{5}(Cntmp)'];
    CRy1 = ones(size(CRx,1),1)+1;
    
    %FOR FA VS CR
    %     [FAx,FAy,CRx,CRy1] = FACRBalance(FAx,CRx);
    %     FCratio(rec,1) = size([FAx],1)/size([FAx; CRx],1);
    %     FCratio(rec,2) = 1-FCratio(rec,1);
    %
    
    switch params
        case 1
            %FOR TOUCH COUNTS
            switch designvars
                case 'counts'
                    hx = [V(rec).touchNum.hit'];
                    FAx = [V(rec).touchNum.FA'];
                    CRx = [V(rec).touchNum.CR'];
                    x = [0:1:10];
                case 'theta'
                    hx = [V(rec).var.hit{1}'];
                    FAx = [V(rec).var.FA{1}'];
                    CRx = [V(rec).var.CR{1}'];
                    x = [min([hx;FAx;CRx])-2:1:max([hx;FAx;CRx])+2];
            end
            
            switch classes
                case 'FAvsCR'
                    firstvar = histc([FAx],x);secondvar = histc([CRx],x);
                case 'lick'
                    firstvar = histc([hx;FAx],x);secondvar = histc([CRx],x);
                case 'gonogo'
                    firstvar = histc([hx],x); secondvar = histc([FAx;CRx],x);
            end
            
            
            figure(12);subplot(2,5,rec)
            h2a=plot(x,firstvar/sum(firstvar),colors{1});h2a.Color(4)=0.5;
            hold on;h2a=plot(x,secondvar/sum(secondvar),colors{2});h2a.Color(4)=0.5;
            
            for db = [1:2]
                ms=cell2mat(Bfeat{rec}.theta);
                coords=mean(reshape(ms(db,:)',2,numIterations),2);
                y= (exp(coords(1)+coords(2)*x)) ./ (1 + exp(coords(1)+coords(2)*x))  ;
                y1=train_predOpt(rec);
                hold on; plot(x,y,['-.o' colors{db}]);
                
                set(gca,'xlim',[min(x) max(x)],'ylim',[0 1]);
                if strcmp(designvars,'counts')
                    xlabel('Counts');
                elseif strcmp(designvars,'theta')
                    xlabel('Theta at Touch');
                end
            end
            if rec ==5
                if strcmp(classes,'gonogo')
                    legend('Go', 'No Go','Go Classifier','No Go Classifier')
                elseif strcmp(classes,'lick')
                    legend('Lick', 'No Lick','Lick Classifier','No Lick Classifier')
                end
            elseif rec == 6
                ylabel('Proportion of Trials / P(TrialType)');
            end
            
            
        case 2
            % %ASIDE: Plotting Decision Boundary for 2 variables
            figure(10);clf;scatter(hx(:,1),hx(:,2),'b')
            hold on;scatter([FAx(:,1);CRx(:,1)],[FAx(:,2);CRx(:,2)],'r')
            xlabel('Theta');ylabel('Velocity')
            
            ms=cell2mat(Bfeat{rec}.theta{1});
            coords=mean(reshape(ms(1,:)',3,10),2);
            
            alldat=filex_whiten([hx;FAx;CRx]);
            figure;scatter(alldat(1:length(hx),1),alldat(1:length(hx),2),'b')
            hold on;scatter(alldat(length(hx):end,1),alldat(length(hx):end,2),'r')
            
            plot_x = [min(alldat(:,1)),  max(alldat(:,1))];
            plot_y = (-1/coords(3)) .* (coords(1) + coords(2).*plot_x - log(train_predOpt{1}(rec)/(1-train_predOpt{1}(rec))));
            hold on;plot(plot_x,plot_y,'k-')
            
        case 3
            % %ASIDE: Plotting Decision Boundary for 3 variables
            
            ms=cell2mat(Bfeat{rec}.theta);
            coords=mean(reshape(ms(1,:)',4,numIterations),2);
            alldat=filex_whiten([hx;FAx;CRx]);
            %             centroid=[mean(alldat(1:length(FAx),:));mean(alldat(length(FAx):end,:))];
            
            figure(10);hold on;subplot(1,2,rec);
            scatter3(alldat(1:length(hx),1),alldat(1:length(hx),2),alldat(1:length(hx),3),'b')
            hold on;scatter3(alldat(length(hx):end,1),alldat(length(hx):end,2),alldat(length(hx):end,3),'r')
            %             hold on;scatter3(centroid(1,1),centroid(1,2),centroid(1,3),'k','linewidth',10)
            %             hold on;scatter3(centroid(2,1),centroid(2,2),centroid(2,3),'b','linewidth',10)
            
            plot_x = [min(alldat(:,1))-2, min(alldat(:,1))-2, max(alldat(:,1))+2, max(alldat(:,1))+2]; %ranges for amplitude
            plot_z = [-3 ,3,3,-3];
            plot_y = (-1/coords(3)) .* (coords(1) + (coords(2).*plot_x) + (coords(4).*plot_z) - log(train_predOpt(rec)/(1-train_predOpt(rec)))); % log p(go trial) to calculate decision boundary
            
            hold on; fill3(plot_x, plot_y, plot_z,'k'); alpha (.4)
            xlabel('Amplitude');ylabel('Midpoint');zlabel('Phase')
            %             legend('Go','Nogo')
            
            ms=cell2mat(Bfeat{rec}.theta);
            optfeat(:,rec)=mean(reshape(ms(1,:)',4,numIterations),2);
            optfeatstd(:,rec)=std(reshape(ms(1,:)',4,numIterations),0,2);
        case 4
            % %ASIDE: Plotting Decision Boundary for 3 variables
            
            
            alldat=filex_whiten([hx;FAx;CRx]);
            centroid=[mean(alldat(1:length(FAx),:));mean(alldat(length(FAx):end,:))];
            
            figure(10);clf;%hold on;subplot(2,3,rec);
            scatter3(alldat(1:length(hx),1),alldat(1:length(hx),2),alldat(1:length(hx),3),'b')
            hold on;scatter3(alldat(length(hx):length([hx;FAx]),1),alldat(length(hx):length([hx;FAx]),2),alldat(length(hx):length([hx;FAx]),3),'g')
            hold on;scatter3(alldat(length([hx;FAx]):end,1),alldat(length([hx;FAx]):end,2),alldat(length([hx;FAx]):end,3),'r')
            
            %             hold on;scatter3(centroid(1,1),centroid(1,2),centroid(1,3),'k','linewidth',10)
            %             hold on;scatter3(centroid(2,1),centroid(2,2),centroid(2,3),'b','linewidth',10)
            for db = [1] %db for go,nogo, hit, fa, cr etc...
                ms=cell2mat(Bfeat{rec}.theta{1});
                coords=mean(reshape(ms(db,:)',4,numIterations),2);
                plot_x = [min(alldat(:,1))-2, min(alldat(:,1))-2, max(alldat(:,1))+2, max(alldat(:,1))+2]; %ranges for amplitude
                plot_z = [-3 ,3,3,-3];
                plot_y = (-1/coords(3)) .* (coords(1) + (coords(2).*plot_x) + coords(4).*plot_z - log(train_predOpt{1}(rec)/(1-train_predOpt{1}(rec)))); % log p(go trial) to calculate decision boundary
                hold on; fill3(plot_x, plot_y, plot_z,'k'); alpha (.4)
                xlabel('Amplitude');ylabel('Midpoint');zlabel('Phase')
            end
            ms=cell2mat(Bfeat{rec}.theta{1});
            optfeat(:,rec)=mean(reshape(ms(1,:)',4,numIterations),2);
            optfeatstd(:,rec)=std(reshape(ms(1,:)',4,numIterations),0,2);
            
            
    end
    
    
end

%% FA:CR DISTRIBUTION
% myC=[0 1 0;1 0 0]; % make a colors list
% H=bar(FCratio, 'stack');
% for k=1:2
%     set(H(k),'facecolor',myC(k,:))
% end
% xlabel('Mouse Number'); ylabel('Proportion of Trials')

