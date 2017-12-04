type = {BV};
for d = 1:length(type) 
U=type{d} ;
clear V F1G RMSEgroup F1Gtree RMSEdecomp 
[V] = classifierWrapper(U);

%% PARAMETERS SETTING

clearvars -except V U BV D SM F1G RMSEgroup RMSEdecomp F1Gtree R PAS POP type d
numIterations = 10;

vars = {'counts'};

for k = 1:length(vars)
designvars = vars{k};
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
    
    [DmatX, DmatY, motorX] = designMatrixBuilder(V(rec),U{rec},designvars,classes,normalization,removal,balance,nanMethod);
    
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
        hold on; plot([0 2],[0 2],'-.k')
        legend('Discrete','Semi Continuous','Continuous','location','northwest')
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


% MAPPING OUT PAS WEIGHTS
if strcmp(designvars,'pas')
    
    pas = abs(optfeat(1:end-1,:));
    normPAS=pas./max(pas);
    if strcmp(U{1}.meta.layer,'SM')
        PAS{2} = normPAS;
    elseif strcmp(U{1}.meta.layer,'BV')
        PAS{3} = normPAS;
    elseif strcmp(U{1}.meta.layer,'D')
        PAS{1} = normPAS;
    end
    
    
%     if ~isempty(PAS{1}) && ~isempty(PAS{2}) &&~isempty(PAS{3})
        PASmean=cell2mat(cellfun(@(x) mean(x,2),PAS,'uniformoutput',0));
        PASstd=cell2mat(cellfun(@(x) std(x,[],2),PAS,'uniformoutput',0));
        figure(90);clf;
        h= filex_barwitherr(PASstd',PASmean');
        ylabel('Normalized Weight')
        set(gca,'xticklabels',{'Discrete','Semi-Continuous','Continuous'},'ytick',[0 .5 1],'ylim',[0 1.2])
        legend('Amplitude','Midpoint','Phase')
        colors = parula(12) ;
        set(gca,'ylim',[0 1.75])
        for d = 1:3
            h(d).FaceColor = [colors(4*d,:)];
        end
%     end
end
% Psychometric Curve Comparison b/t Model and Mouse
if strcmp(classes,'lick')
    [RMSE] = rebuildPsycho(U,V,train_motorPlick);
    suptitle([U{rec}.meta.layer ' ' designvars ' ' classes])
    print(figure(5),'-dtiff',['Z:\Users\Jon\Projects\Characterization\' U{rec}.meta.layer '\Figures\'  U{rec}.meta.layer '_' classes '_' designvars '_LOGPSYCHO'])
    
    if strcmp(designvars,'counts')
        RMSEgroup(:,1) = RMSE;
        
    elseif strcmp(designvars,'theta')
        RMSEgroup(:,2) = RMSE;
        
    elseif strcmp(designvars,'ubered')
        RMSEgroup(:,3) = RMSE;
        
    end
end
end
end