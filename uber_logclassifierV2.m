


clear V
clear F1G RMSEgroup F1Gtree RMSEdecomp 
[V] = classifierWrapper(U);

%% PARAMETERS SETTING

clearvars -except V U BV D SM F1G RMSEgroup RMSEdecomp F1Gtree R PAS
numIterations = 10;

designvars = 'ubered';
% 1) 'theta' 2) 'pas' (phase amp midpoint) 3) 'counts' $) 'ubered'
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

%% LOG CLASSIFIER
clear Acc
accprop=cell(1,length(V));
for rec = 1:length(V)
    
    [DmatX, DmatY, motorX] = designMatrixBuilder(V(rec),U{rec},designvars,classes,normalization,removal,balance);
    
    g1 = DmatY(DmatY == 1);
    g2 = DmatY(DmatY == 2);
    
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
                normPAS = [1:length(tmpDmatY)]';
                normPAS(train)=[];
                
                [thetas,cost,~] = ML_oneVsAll(tmpDmatX(train,:),tmpDmatY(train,:),numel(unique(DmatY)),0);
                Bfeat{rec}.theta{reit}=thetas;
                
                [pred,opt_thresh(reit),prob]=ML_predictOneVsAll(thetas,tmpDmatX(normPAS,:)...
                    ,tmpDmatY(normPAS,:),'Max');
                Acc(reit)= mean(double(pred == tmpDmatY(normPAS))) * 100;
                F1s(reit,:) = F1score(pred,tmpDmatY(normPAS),2);
                
                if strcmp(classes,'lick')
                    motorPlick= [motorPlick;tmpMotorX(normPAS) prob];
                end
                accprop{rec}=[accprop{rec} ; pred tmpDmatY(normPAS)];
                
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
        
    end
    if strcmp(classes,'lick')
        train_motorPlick{rec} = motorPlick;
    end
    train_F1s(rec,:) = nanmean(nansum(F1s,2));
    trainF1sstd(rec,:)=nanstd(nansum(F1s,2)); 
    train_predOpt(rec)=mean(opt_thresh); %used for dboundaries
    train_Acc(rec) = mean(Acc); 
    train_std(rec) = std(Acc);


end


%Plotting F1 score for Log Classifier
if strcmp(designvars,'theta')
    F1G{1} = train_F1s ;
elseif strcmp(designvars,'counts')
    F1G{2} = train_F1s ;
    figure(20);hold on;f=scatter(F1G{1},F1G{2},'filled');
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

% Plotting Model accuracy for Log Classifier
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
        tmp=find(guesses(:,2)==i); %find predictions of class i 
        predprop(recs,i)=mean(guesses(tmp,1)==i)*100;
    end
end
hold on; scatter(1:length(V),predprop(:,1),'b');
scatter(1:length(V),predprop(:,2),'r')
legend('Model Accuracy','Lick Prediction Accuracy','No Lick Prediction Accuracy','location','southeast')
title(['Log' U{rec}.meta.layer ' ' designvars ' ' classes])
print(figure(22),'-dtiff',['Z:\Users\Jon\Projects\Characterization\' U{rec}.meta.layer '\Figures\'  U{rec}.meta.layer '_' classes '_' designvars '_0touch_removal_' removal])

%Plotting Weights for Log Classifier 
if strcmp(designvars,'ubered') || strcmp(designvars,'pas')
for rec = 1:length(V)
    ms=cell2mat(Bfeat{rec}.theta);
    optfeat(:,rec)=mean(reshape(ms(1,:)',size(DmatX,2)+1,numIterations),2);
    optfeatstd(:,rec)=std(reshape(ms(1,:)',size(DmatX,2)+1,numIterations),0,2);
end
optfeatstd = [optfeatstd(2:end,:);optfeatstd(1,:)];
optfeat = [optfeat(2:end,:);optfeat(1,:)];
figure(92);h=filex_barwitherr(optfeatstd',optfeat');

ylabel('Relative Weight');
if strcmp(designvars,'pas')
    h(4).FaceColor = [.8 .8 .8]; %make bias weight color gray
    legend('Amplitude at Touch','Midpoint at Touch','Phase at Touch','Bias','Location','northeast')
elseif strcmp(designvars,'ubered')
    h(6).FaceColor = [.8 .8 .8]; %make bias weight color gray
    legend('Theta at Touch','Touch Count','Previous Trial Lick','2 Previous T Lick', '3 Previous T Lick','Bias','Location','northeast')
end

print(figure(92),'-dtiff',['Z:\Users\Jon\Projects\Characterization\' U{rec}.meta.layer '\Figures\'  U{rec}.meta.layer '_' classes '_' designvars '_WEIGHTS'])
end


%% MAPPING OUT PAS WEIGHTS
pas = abs(optfeat(1:end-1,:));
normPAS=pas./max(pas);
if strcmp(U{1}.meta.layer,'SM')
    PAS{2} = normPAS;
elseif strcmp(U{1}.meta.layer,'BV')
    PAS{3} = normPAS;
 elseif strcmp(U{1}.meta.layer,'D')
    PAS{1} = normPAS;
end


if ~isempty(PAS{1}) && ~isempty(PAS{2}) &&~isempty(PAS{3})
    PASmean=cell2mat(cellfun(@(x) mean(x,2),PAS,'uniformoutput',0));
    PASstd=cell2mat(cellfun(@(x) std(x,[],2),PAS,'uniformoutput',0));
    figure(90);clf;
    h= filex_barwitherr(PASstd',PASmean');
    ylabel('Normalized Weight')
    set(gca,'xticklabels',{'Discrete','Semi-Continuous','Continuous'},'ytick',[0 .5 1],'ylim',[0 1.2])
    legend('Amplitude','Midpoint','Phase')
    colors = parula(12) ;
    for d = 1:3
        h(d).FaceColor = [colors(4*d,:)];
    end
end

% Psychometric Curve Comparison b/t Model and Mouse
if strcmp(classes,'lick')
[RMSE] = rebuildPsycho(U,V,train_motorPlick);
suptitle([U{rec}.meta.layer ' ' designvars ' ' classes])
print(figure(5),'-dtiff',['Z:\Users\Jon\Projects\Characterization\' U{rec}.meta.layer '\Figures\'  U{rec}.meta.layer '_' classes '_' designvars '_LOGPSYCHO'])

    if strcmp(designvars,'counts')
        RMSEgroup(rec,1) = RMSE;

    elseif strcmp(designvars,'theta')
        RMSEgroup(rec,2) = RMSE;

    elseif strcmp(designvars,'ubered')
        RMSEgroup(rec,3) = RMSE;

    end
    
    if size(RMSEgroup,2)==3

        if strcmp(U{1}.meta.layer,'BV')
            R{3} = RMSEgroup;
        elseif strcmp(U{1}.meta.layer,'SM')
            R{2} = RMSEgroup;
        elseif strcmp(U{1}.meta.layer,'D')
            R{1} = RMSEgroup;
        end
    end
end
%% Visualization of the decision boundaries.
params =1;
[V] = classifierWrapper(U);
colors = {'b','r'}

for rec = 1:length(U)
 
    
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

%% TREE BAG CLASSIFIER 
reitAcc = cell(1,length(U));
reitPred = cell(1,length(U));
reitReal = cell(1,length(U));
reitMotor = cell(1,length(U));
reitBFeat = cell(1,length(U));
groupF1 = nan(length(U),2); 
for rec = 1:length(U)
    
   for reit = 1:numIterations
    [DmatX, DmatY , motorX] = designMatrixBuilder(V(rec),U{rec},designvars,classes,normalization,removal,balance,nanMethod);
    
    [bagPred ,bagReal, motorPos, treeMdl] = uber_baggingClassifier(DmatX,DmatY,motorX);
    
    
    reitAcc{rec}(reit) = mean(bagPred(:,1)==bagReal);
    reitPred{rec}=[reitPred{rec} ; bagPred(:,1)];
    reitReal{rec}=[reitReal{rec} ; bagReal];
    reitMotor{rec} = [reitMotor{rec} ; motorPos bagPred(:,2)];
    reitBFeat{rec} = [reitBFeat{rec} ; treeMdl.OOBPermutedVarDeltaError];
    
    
   end
   groupF1(rec,:) = F1score(reitPred{rec},reitReal{rec},2);
end

%Plotting average test accuracy 
figure(36);clf
errorbar(1:length(U),cellfun(@mean,reitAcc)*100,cellfun(@std,reitAcc)*100,'ko')
xlabel('Mouse Number');ylabel('% Accuracy of Model')
set(gca,'ylim',[0 100])
xlim([.5 length(V)+.5])

%plotting TRUE POSITIVES (GUESSING class 1 and then it actually being class
%1)
predprop = zeros(length(V),3);
for recs = 1:length(V)
    guesses=reitPred{recs};
    real = reitReal{recs};
    for i = 1:3
        tmp=find(guesses==i);
        predprop(recs,i)=mean(real(tmp)==i)*100;
    end
end
hold on; scatter(1:length(V),predprop(:,1),'b');
scatter(1:length(V),predprop(:,2),'r')
legend('Model Accuracy','Lick Prediction Accuracy','No Lick Prediction Accuracy','location','southeast')
title(['TreeBag ' U{rec}.meta.layer ' ' designvars ' ' classes])
 print(figure(36),'-dtiff',['Z:\Users\Jon\Projects\Characterization\' U{rec}.meta.layer '\Figures\'  U{rec}.meta.layer '_' classes '_' designvars '_0touch_removal_' removal '_TREEBAG'])

%Plotting Feature Importance for Trees
if strcmp(designvars,'ubered') || strcmp(designvars,'pas')
    treeFeat = cell2mat(cellfun(@mean,reitBFeat,'UniformOutput',0)');
    treeStd = cell2mat(cellfun(@std,reitBFeat,'UniformOutput',0)');
    figure(323);filex_barwitherr(treeStd,treeFeat)
    ylabel('Tree Feature Importance');
    if strcmp(designvars,'pas')
        legend('Amplitude at Touch','Midpoint at Touch','Phase at Touch','Location','northeast')
    elseif strcmp(designvars,'ubered')
        legend('Theta at Touch','Touch Count','Previous Trial Lick','2 Previous T Lick', '3 Previous T Lick','Location','northeast')
    end
    print(figure(323),'-dtiff',['Z:\Users\Jon\Projects\Characterization\' U{rec}.meta.layer '\Figures\'  U{rec}.meta.layer '_' classes '_' designvars '_0touch_removal_' removal '_TREEBAG_WEIGHTS'])

end

%Plotting F1 score for Trees
for d = 1:length(V)
    F1group(d) = nansum(F1score(reitPred{d},reitReal{d},2));

end

if strcmp(designvars,'theta')
    F1Gtree{1} = F1group;
elseif strcmp(designvars,'counts')
    F1Gtree{2} = F1group;
        figure(32);hold on;f=scatter(F1Gtree{1},F1Gtree{2},'filled');
    if strcmp(U{rec}.meta.layer,'D')
        f.CData = [rgb('DarkGreen')];
    elseif strcmp(U{rec}.meta.layer,'SM')
        f.CData = [rgb('DarkMagenta')];
    elseif strcmp(U{rec}.meta.layer,'BV')
        f.CData = [rgb('DarkTurquoise')];
    end
    set(gca,'xlim',[0 2],'ylim',[0 2],'ytick',0:.5:2);
    xlabel('F1 Score Theta');ylabel('F1 Score Touch Count');
        hold on; plot([0 2],[0 2],'-.k')
        legend('Discrete','Semi Continuous','Continuous','location','northwest')
 
end

% Psychometric Curve Comparison b/t Model and Mouse
[RMSE] = rebuildPsycho(U,V,reitMotor);
suptitle([U{rec}.meta.layer ' ' designvars ' ' classes])
print(figure(5),'-dtiff',['Z:\Users\Jon\Projects\Characterization\' U{rec}.meta.layer '\Figures\'  U{rec}.meta.layer '_' classes '_' designvars '_TREEPSYCHO'])

    if strcmp(designvars,'counts')
        RMSEgroup(rec,1) = RMSE;

    elseif strcmp(designvars,'theta')
        RMSEgroup(rec,2) = RMSE;

    elseif strcmp(designvars,'ubered')
        RMSEgroup(rec,3) = RMSE;

    end
    
    if size(RMSEgroup,2)==3

        if strcmp(U{1}.meta.layer,'BV')
            R{3} = RMSEgroup;
        elseif strcmp(U{1}.meta.layer,'SM')
            R{2} = RMSEgroup;
        elseif strcmp(U{1}.meta.layer,'D')
            R{1} = RMSEgroup;
        end
    end
