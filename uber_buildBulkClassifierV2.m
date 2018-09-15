clearvars -except SM BV
type = {SM,BV};
colors = {'Gold','DarkTurquoise'};
PAS = cell(1,length(type));
for d =2
    U=type{d} ;
    clear V F1G RMSEgroup F1Gtree RMSEdecomp
    [V] = classifierWrapper(U);
    
    %% PARAMETERS SETTING
    
%     clearvars -except V U BV D SM POPv3 type d  countsF1 thetaF1 uberedF1 timingF1 motorF1 uberedfeat countsfeat thetafeat motorfeat RMSEpop timingfeat pasfeat dtimefeat
    numIterations = 20;
    
    
    %     vars = {'counts','theta','timing','motor','ubered'};
    vars = {'theta','pas','decompTime'};
    for k = 1:3
        designvars = vars{k};
        % 1) 'theta' 2) 'pas' (phase amp midpoint) 3) 'counts' 4) 'ubered'
        % 5) 'timing' 6) 'motor' 7) 'decompTime'
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
        %% LOG CLASSIFIER
        clear Acc Bfeat
        accprop=cell(1,length(V));
        for rec = 1:length(V)
            
            %     [DmatX, DmatY, motorX] = designMatrixBuilder(V(rec),U{rec},designvars,classes,normalization,removal,balance,nanMethod);
            % V2 USES THETA COLUMN FROM UBERED MODEL
            [DmatX, DmatY, motorX] = designMatrixBuilderv2(V(rec),U{rec},designvars,classes,normalization,removal,balance,nanMethod,dropempty);
            
            %             DmatX = motorX;
            
            g1 = DmatY(DmatY == 1);
            g2 = DmatY(DmatY == 2);
            
            clear opt_thresh
            motorPlick = [];
            txp = [];
            for reit = 1:numIterations
                rando = randperm(length(DmatX));
                tmpDmatX=DmatX(rando,:);tmpDmatY=DmatY(rando,:);
                %                 if strcmp(classes,'lick')
                tmpMotorX = motorX(rando,:);
                %                 end
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
                        
                        %                         if strcmp(classes,'lick')
                        motorPlick= [motorPlick;tmpMotorX(normPAS) prob];
                        %                         end
                        accprop{rec}=[accprop{rec} ; pred tmpDmatY(normPAS)];
                        
                        txp = [txp ; tmpDmatY(normPAS) pred];
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
            
            poptxp{rec} = txp;
            popDmatX{rec} = DmatX; 
            train_motorPlick{rec} = motorPlick; %used for psycho curves
            train_predOpt(rec)=mean(opt_thresh); %used for dboundaries
            train_Acc(rec) = mean(Acc);
            train_std(rec) = std(Acc);
            
            
            
        end
        
        
        %Raw truths and predictions: use this to calculate MCC, F1 scores,
        %etc...
        if strcmp(designvars,'theta')
            thetatxp{d}=poptxp;
        elseif strcmp(designvars,'counts')
            countstxp{d}=poptxp;
        elseif strcmp(designvars,'timing')
            timingtxp{d}=poptxp;
        elseif strcmp(designvars,'ubered')
            uberedtxp{d}=poptxp;
        elseif strcmp(designvars,'motor')
            motortxp{d} =poptxp;
        elseif strcmp(designvars,'pas')
            pastxp{d} =poptxp;
        elseif strcmp(designvars,'decompTime')
            dtimetxp{d} =poptxp;
        end
        
        
        %
        if strcmp(designvars,'theta')
            thetaDmatX{d}=popDmatX;
        elseif strcmp(designvars,'counts')
            countsDmatX{d}=popDmatX;
        elseif strcmp(designvars,'timing')
            timingDmatX{d}=popDmatX;
        elseif strcmp(designvars,'ubered')
            uberedDmatX{d}=popDmatX;
        elseif strcmp(designvars,'motor')
            motorDmatX{d} =popDmatX;
        elseif strcmp(designvars,'pas')
            pasDmatX{d} =popDmatX;
        elseif strcmp(designvars,'decompTime')
            dtimeDmatX{d} =popDmatX;
        end
        
        
        % Plotting Weights for Log Classifier
        if strcmp(designvars,'ubered') || strcmp(designvars,'pas') || strcmp(designvars,'decompTime')
            optfeat = zeros(size(DmatX,2)+1,length(V));
            optfeatstd = zeros(size(DmatX,2)+1,length(V));
            for rec = 1:length(V)
                ms=cell2mat(Bfeat{rec}.theta);
                optfeat(:,rec)=mean(reshape(ms(1,:)',size(DmatX,2)+1,numIterations),2);
                optfeatstd(:,rec)=std(reshape(ms(1,:)',size(DmatX,2)+1,numIterations),0,2);
            end
            optfeatstd = [optfeatstd(2:end,:);optfeatstd(1,:)];
            optfeat = [optfeat(2:end,:);optfeat(1,:)];

            if strcmp(designvars,'ubered')
                uberedfeat{d} = optfeat;
            elseif strcmp(designvars,'pas')
                pasfeat{d} = optfeat;
            elseif strcmp(designvars,'decompTime')
                dtimefeat{d} = optfeat;
            end
            
            %             figure(92);h=filex_barwitherr(optfeatstd',optfeat');

        else
            for rec = 1:length(V)
                ms=cell2mat(Bfeat{rec}.theta);
                optfeat(:,rec)=mean(reshape(ms(1,:)',size(DmatX,2)+1,numIterations),2);
                optfeatstd(:,rec)=std(reshape(ms(1,:)',size(DmatX,2)+1,numIterations),0,2);
            end
            
            optfeatstd = [optfeatstd(2:end,:);optfeatstd(1,:)];
            optfeat = [optfeat(2:end,:);optfeat(1,:)];
            if strcmp(designvars,'theta')
                thetafeat{d} = optfeat;
            elseif strcmp(designvars,'counts')
                countsfeat{d} = optfeat;
            elseif strcmp(designvars,'timing')
                timingfeat{d} = optfeat;
            elseif strcmp(designvars,'motor')
                motorfeat{d} = optfeat;
            end
            
        end
        
        
        
        %         % MAPPING OUT PAS WEIGHTS
        %         if strcmp(designvars,'pas')
        %
        %             pas = abs(optfeat(1:end-1,:));
        %             normPAS=pas./max(pas);
        %             if strcmp(U{1}.meta.layer,'SM')
        %                 PAS{1} = normPAS;
        %             elseif strcmp(U{1}.meta.layer,'BV')
        %                 PAS{2} = normPAS;
        %             elseif strcmp(U{1}.meta.layer,'D')
        %                 %         PAS{1} = normPAS;
        %             end
        %
        %
        %             %     if ~isempty(PAS{1}) && ~isempty(PAS{2}) &&~isempty(PAS{3})
        %             if ~isempty(PAS{1}) && ~isempty(PAS{2})
        %                 PASmean=cell2mat(cellfun(@(x) mean(x,2),PAS,'uniformoutput',0));
        %                 PASstd=cell2mat(cellfun(@(x) std(x,[],2),PAS,'uniformoutput',0));
        %                 figure(90);clf;
        %                 h= filex_barwitherr(PASstd',PASmean');
        %                 ylabel('Normalized Weight')
        %                 set(gca,'xticklabels',{'Semi-Continuous','Continuous'},'ytick',[0 .5 1],'ylim',[0 1.2])
        %                 legend('Amplitude','Midpoint','Phase')
        %                 colors = parula(12) ;
        %                 set(gca,'ylim',[0 1.75])
        %                 for d = 1:2
        %                     h(d).FaceColor = [colors(4*d,:)];
        %                 end
        %             end
        %         end
        
        
        %         Psychometric Curve Comparison b/t Model and Mouse
        %         if strcmp(classes,'lick')
        %             [RMSE] = rebuildPsychoflip(U,V,train_motorPlick);
        %             suptitle([U{rec}.meta.layer ' ' designvars ' ' classes])
        % %             print(figure(5),'-depsc',['C:\Users\jacheung\Dropbox\LocalizationBehavior\Figures\Parts\'  U{rec}.meta.layer '_' classes '_' designvars '_LOGPSYCHO'])
        %
        %             if strcmp(designvars,'counts')
        %                 RMSEgroup(:,1) = RMSE;
        %             elseif strcmp(designvars,'theta')
        %                 RMSEgroup(:,2) = RMSE;
        %             elseif strcmp(designvars,'timing')
        %                 RMSEgroup(:,3) = RMSE;
        %             elseif strcmp(designvars,'ubered')
        %                 RMSEgroup(:,4) = RMSE;
        %             end
        %             RMSEpop{d} = RMSEgroup;
        %         end
        
    end
    %% Visualization of the decision boundaries.
    params =3;
    colors = {'b','r'};
    figure(12);clf
    for rec = 1:length(U)
        switch params
            case 1
    
                [DmatX, DmatY, motorX] = designMatrixBuilderv2(V(rec),U{rec},designvars,classes,normalization,removal,balance,nanMethod,dropempty);
    %             x=linspace(min(DmatX),max(DmatX),max(DmatX)-min(DmatX)+1);
    %             x=linspace(min(DmatX),max(DmatX),11);
                x=linspace(0,750,25);
                firstvar = histc(DmatX(DmatY==1),x);
                secondvar = histc(DmatX(DmatY==2),x);
    
    
    %             switch classes
    %                 case 'lick'
    %                     firstvar = histc([hx;FAx],x);secondvar = histc([CRx],x);
    %                 case 'gonogo'
    %                     firstvar = histc(DmatX(DmatY==1),x); secondvar = histc(DmatX(DmatY==2),x);
    %             end
    %
    
                figure(12);subplot(2,5,rec)
                bar(x,firstvar/sum(firstvar),colors{1});
                hold on;bar(x,secondvar/sum(secondvar),colors{2});
    
                for db = [1:2]
                        ms=cell2mat(Bfeat{rec}.theta);
                    coords=mean(reshape(ms(db,:)',2,numIterations),2);
                    y= (exp(coords(1)+coords(2)*x)) ./ (1 + exp(coords(1)+coords(2)*x))  ;
                    y1=train_predOpt(rec);
                    hold on; plot(x,y,['-.' colors{db}]);
    
                    set(gca,'xlim',[min(x) max(x)],'ylim',[0 1]);
                    if strcmp(designvars,'counts')
                        xlabel('Counts');
                    elseif strcmp(designvars,'theta')
                        xlabel('Theta at Touch');
                    end
                end
                set(gca,'ytick',0:.5:1,'xtick',0:250:750)
    
    
    %             if rec ==5
    %                 if strcmp(classes,'gonogo')
    %                     legend('Go', 'No Go','Go Classifier','No Go Classifier')
    %                 elseif strcmp(classes,'lick')
    %                     legend('Lick', 'No Lick','Lick Classifier','No Lick Classifier')
    %                 end
    %             elseif rec == 6
    %                 ylabel('Proportion of Trials / P(TrialType)');
    %             end
       case 3
            % %ASIDE: Plotting Decision Boundary for 3 variables
            


            [DmatX, DmatY, motorX] = designMatrixBuilderv2(V(rec),U{rec},designvars,classes,normalization,removal,balance,nanMethod,dropempty);
%             coords = optfeat(:,rec);
            ms=cell2mat(Bfeat{rec}.theta);
            coords=mean(reshape(ms(1,:)',4,numIterations),2);
            
            figure(10);clf
            scatter3(DmatX(DmatY==1,1),DmatX(DmatY==1,2),DmatX(DmatY==1,3),'m.')
            hold on;scatter3(DmatX(DmatY==2,1),DmatX(DmatY==2,2),DmatX(DmatY==2,3),'k.')
            
            %             hold on;scatter3(centroid(1,1),centroid(1,2),centroid(1,3),'k','linewidth',10)
            %             hold on;scatter3(centroid(2,1),centroid(2,2),centroid(2,3),'b','linewidth',10)
            
            plot_x = [min(DmatX(:,1))-2, min(DmatX(:,1))-2, max(DmatX(:,1))+2, max(DmatX(:,1))+2]; %ranges for amplitude
            plot_z = [-3 ,3,3,-3];
            plot_y = (-1/coords(3)) .* (coords(1) + (coords(2).*plot_x) + (coords(4).*plot_z) - log(train_predOpt(rec)/(1-train_predOpt(rec)))); % log p(go trial) to calculate decision boundary
            
             hold on; fill3(plot_x, plot_y, plot_z,'k');
            xlabel('Amplitude');ylabel('Midpoint');zlabel('Phase')
            set(gca,'xlim',[-3 3],'ylim',[-3 3],'zlim',[-3 3])
            %             legend('Go','Nogo')
            

        end
    
    
    end
end


