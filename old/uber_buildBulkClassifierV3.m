clearvars -except SM BV
type = {SM,BV};
colors = {'Gold','DarkTurquoise'};
PAS = cell(1,length(type));
for d =2
    U=type{d} ;
    clear V F1G RMSEgroup F1Gtree RMSEdecomp
    [V] = classifierWrapper(U);
    
    
    %% PARAMETERS SETTING
    
    cvKfold = 5;

    %     vars = {'midpoint','amp','phase','mpamp','mpphase','ampphase'};
    %     vars = {'midpoint','amp','phase','mpamp','mpphase','ampphase','countsmidpoint','countsamp','countsphase','countsmpamp','countsmpphase','countsampphase'};
    
%     vars = {'countsBinary','counts'}; %Fig 5
    %     vars = {'kappa','timing','whiskOnset','counts','radialD','theta'}; %Fig 7
%     vars = {'motor','kappa','timing','whiskOnset','counts','radialD','theta','uberedRadial'}; %Fig 9
%     vars = {'theta','pas','decompTime'} %Fig 10
%     vars = {'timeTotouch','onsetTheta','velocity','phase','amp','midpoint','theta'} %Fig11
%     vars = {'timeTotouch','onsetTheta','velocity','decompTime'};
    vars = {'uberedRadial'}

    
    for k = 1:length(vars)
        designvars = vars{k};
        % 1) 'theta' 2) 'pas' (phase amp midpoint) 3) 'counts' 4) 'ubered'
        % 5) 'timing' 6) 'motor' 7) 'decompTime' OR ,'kappa'
        % 'timeTotouch','onsetTheta','velocity','Ivelocity' OR 'phase','amp','midpoint'
        classes = 'lick';
        % 1) 'gonogo' 2) 'FAvsCR' 3) 'lick'
        
        % Only for 'ubered' or 'pas'
        normalization = 'meanNorm';
        % 1) 'whiten'  2) 'meanNorm' 3)'none'
        
        %Regularization method,
        regMethod = 'lasso';
        lambda = loadLambda;
        % 1) 'lasso' or L1 regularization;
        % 2) 'ridge' or L2 regularziation;
        
        biasClose = 'no';
        distance_round_pole =2; %inmm
        
        
        nanMethod = 'random';
        % doesnt really matter because we end up dropping non touch trials
        % anyways
        % 1) random (resample NaN vars from all touches)
        % 2) peakPro (replace NaN using var from max protraction)
        % 3) resampSameD (replace NaN using vas from touches in same trial type)
        
        dropempty = 'yes';
        % 1) 'yes' = drop trials with 0 touches
        % 2) 'no' = keep all trials
        %% LOG CLASSIFIER
        clear indivMedians
        for z = 1:length(lambda)
            clear Bfeat
            accprop=cell(1,length(V));
            for rec = 1:10
                
                [DmatX, DmatY, motorX] = designMatrixBuilderv3(V(rec),U{rec},designvars,classes,normalization,nanMethod,dropempty);
               
                
                if strcmp(biasClose,'yes')
                    mean_norm_motor = motorX - mean(BV{rec}.meta.ranges);
                    close_trials = find(abs(mean_norm_motor)<distance_round_pole*10000);
                    DmatX = DmatX(close_trials,:);
                    DmatY = DmatY(close_trials,:);
                end
                
                g1 = DmatY(DmatY == 1);
                g2 = DmatY(DmatY == 2);
                
                g1cvInd = crossvalind('kfold',length(g1),5);
                g2cvInd = crossvalind('kfold',length(g2),5);
                
                % shuffle permutations of cv indices
                g1cvInd = g1cvInd(randperm(length(g1cvInd)));
                g2cvInd = g2cvInd(randperm(length(g2cvInd)));
                
                selInds = [g1cvInd ;g2cvInd];
                
                clear opt_thresh
                motorPlick = [];
                txp = [];
                cost = nan(cvKfold,1);
                for u=1:cvKfold
                    testY = [DmatY(selInds==u)];
                    testX = [DmatX(selInds==u,:)];
                    trainY = [DmatY(~(selInds==u))];
                    trainX = [DmatX(~(selInds==u),:)];
                    
                    if size(DmatX,2)<2 %prediction using single features arent lasso'd
                        lambda = 0 ;
                    end
                    [thetas,cost,~] = ML_oneVsAll(trainX, trainY, numel(unique(DmatY)), lambda(z), regMethod);
                    Bfeat{rec}.theta{u}=thetas;
                   
                    [pred,opt_thresh(u),prob]=ML_predictOneVsAll(thetas,testX,testY,'Max');
                    
                    motorPlick= [motorPlick;motorX(selInds==u) prob];
                    txp = [txp ; testY pred];
                end
                
                poptxp{rec} = txp;
                popDmatX{rec} = DmatX;
                train_motorPlick{rec} = motorPlick; %used for psycho curves
                train_predOpt(rec)=mean(opt_thresh); %used for dboundaries
            end
            
            %Raw truths and predictions: use this to calculate MCC, F1 scores,
            %etc...
            trueXpreds{d}.(vars{k}) = poptxp;
            lambdatxp{z} = poptxp;
            
            
            clear feats
            clear oddsR
            if strcmp(designvars,'ubered') || strcmp(designvars,'uberedRadial')|| strcmp(designvars,'pas') || strcmp(designvars,'decompTime')
                
                for rec = 1:length(Bfeat)
                    ms=cell2mat(Bfeat{rec}.theta);
                    feats{rec} = reshape(ms(1,:)',size(DmatX,2)+1,cvKfold);
                end
                
                for g = 1:length(Bfeat)
                    for b = 1:cvKfold
                        xv = feats{g}(2:end,b);
                        pn = xv./abs(xv);
                        
                        pos = exp(xv).*(pn==1);
                        neg = exp(-xv).*(pn==-1);
                        pos(pos==0)=neg(~(neg==0));
                        
                        oddsR{g}(:,b) = pos.*pn;
                    end
                    
                    oddsR{g}(isinf(oddsR{g})) = nan;
                    meanOdds(g,:) = nanmean(oddsR{g},2);
                end
                
                clear normoddsRmeans
                for g=1:length(oddsR)
                    normoddsRmeans(:,g) = nanmedian(normalize_var(abs(oddsR{g}),0,1) .* (oddsR{g}./abs(oddsR{g})) ,2);
                    oddsRmedians(:,g) = nanmedian(abs(oddsR{g}) .* (oddsR{g}./abs(oddsR{g})),2);
                    indivMedians{g}(z,:) = nanmedian(abs(oddsR{g}) .* (oddsR{g}./abs(oddsR{g})),2);
                end
                

%                 oddsRmeans = cell2mat(cellfun(@(x) median(x,2),oddsR,'uniformoutput',0));
%                 oddsRstd = cell2mat(cellfun(@(x) std(x,[],2),oddsR,'uniformoutput',0));
%                 normoddsRmeans = normalize_var(abs(oddsRmeans),0,1).*(oddsRmeans./abs(oddsRmeans));
                
                
                figure(488);clf
                for i = 1:length(Bfeat)
                    %subplot(5,2,k);
                    %                 figure(499);
                    %                 hold on; plot(1:size(oddsRmeans,1),oddsRmeans(:,k),'color',[.8 .8 .8])
                    figure(488);
                    hold on; plot(1:size(normoddsRmeans,1),normoddsRmeans(:,i),'color',[.8 .8 .8])
                    hold on;plot([1 size(DmatX,2)],[0 0],'-.k')
                    set(gca,'xlim',[1 size(DmatX,2)])
                end


                figure(488);
                hold on; plot(1:size(DmatX,2),mean(normoddsRmeans,2),'-ko')
                            set(gca,'xtick',1:size(DmatX,2),'xticklabel',{'kappa','cueTime','onsetTime','counts','radial','angle'},'ytick',[-1:.5:1],'ylim',[-1 1])
                %             set(gca,'xtick',1:size(DmatX,2),'xticklabel',{'amp','mp','phase'},'ytick',[-1:.5:1],'ylim',[-1 1])
%                 set(gca,'xtick',1:size(DmatX,2),'xticklabel',{'time','startT','vel'},'ytick',[-1:.5:1],'ylim',[-1 1])
                
                
                
                title('relative odds')
                ylabel('odds lick per one std increase')
                
                
                [p,~,stats]=anova1(normoddsRmeans',[],'off');
%                 multcompare(stats,[],'off')
                
                keepraw(z,:) = nanmedian(oddsRmedians,2);
                keep(z,:) = mean(normoddsRmeans,2);
                keepstd(z,:) = std(normoddsRmeans,[],2);
                lambdatxp{z} = poptxp;
            end
%             
            
            
        end
        %                 hold on; plot(1:4,mean(oddsRmeans,2)-std(oddsRmeans,[],2),'-ro')
        %                  hold on; plot(1:4,mean(oddsRmeans,2)+std(oddsRmeans,[],2),'-ro')
        %                 set(gca,'yscale','log')
        
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
        
        
        %                 Psychometric Curve Comparison b/t Model and Mouse
        %         if strcmp(classes,'lick')
        %             rebuildPsychoflip_v2(U,V,train_motorPlick);
        %             suptitle([U{rec}.meta.layer ' ' designvars ' ' classes])
        % %             print(figure(5),'-depsc',['C:\Users\jacheung\Dropbox\LocalizationBehavior\Figures\Parts\'  U{rec}.meta.layer '_' classes '_' designvars '_LOGPSYCHO'])
        %
        %         end
        
    end
    %% Visualization of the decision boundaries.
    colors = {'b','r'};
    figure(12);clf
    for rec = 3
       [DmatX, DmatY, motorX] = designMatrixBuilderv3(V(rec),U{rec},designvars,classes,normalization,nanMethod,dropempty);   %             x=linspace(min(DmatX),max(DmatX),max(DmatX)-min(DmatX)+1);
        params = size(DmatX,2);
        switch params
            case 1
                
               x=linspace(min(DmatX),max(DmatX),11);
                
                %                DmatX = DmatX./33;
%                 x=linspace(0,60,15);
                %                 x= linspace(-.05,.05,51);
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
                    coords=mean(reshape(ms(db,:)',2,cvKfold),2);
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
                %                  set(gca,'ytick',0:.5:1,'xtick',0:250:750)
                set(gca,'xlim',[min(DmatX) max(DmatX)])
                
                
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
                
                
                
                [DmatX, DmatY, motorX] = designMatrixBuilderv3(V(rec),U{rec},designvars,classes,normalization,nanMethod,dropempty);%             coords = optfeat(:,rec);
                ms=cell2mat(Bfeat{rec}.theta);
                coords=mean(reshape(ms(1,:)',4,cvKfold),2);
                
                figure(10);clf
                scatter3(DmatX(DmatY==1,1),DmatX(DmatY==1,2),DmatX(DmatY==1,3),'m.')
                hold on;scatter3(DmatX(DmatY==2,1),DmatX(DmatY==2,2),DmatX(DmatY==2,3),'k.')
                
                %             hold on;scatter3(centroid(1,1),centroid(1,2),centroid(1,3),'k','linewidth',10)
                %             hold on;scatter3(centroid(2,1),centroid(2,2),centroid(2,3),'b','linewidth',10)
                
                plot_x = [min(DmatX(:,1))-2, min(DmatX(:,1))-2, max(DmatX(:,1))+2, max(DmatX(:,1))+2]; %ranges for amplitude
                plot_z = [-3 ,3,3,-3];
                plot_y = (-1/coords(3)) .* (coords(1) + (coords(2).*plot_x) + (coords(4).*plot_z) - log(train_predOpt(rec)/(1-train_predOpt(rec)))); % log p(go trial) to calculate decision boundary
                
                hold on; fill3(plot_x, plot_y, plot_z,'k');
                
                set(gca,'xlim',[min(DmatX(:,1)) max(DmatX(:,1))],'ylim',[min(DmatX(:,2)) max(DmatX(:,2))],'zlim',[min(DmatX(:,3)) max(DmatX(:,3))])
                if strcmp(designvars,'pas')
                    xlabel('Amplitude');ylabel('Midpoint');zlabel('Phase')
                elseif strcmp(designvars,'decompTime')
                    xlabel('time to touch');ylabel('whisk angle onset');zlabel('velocity')
                end
                %             legend('Go','Nogo')
                
                
        end
        
        
    end
end


