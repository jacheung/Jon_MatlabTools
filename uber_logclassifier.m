clear; close all;clc;

load('Z:\Users\Jon\Projects\Characterization\BV')
U=BV;
clear V
clear F1G
[V] = classifierWrapper(U);
%% PARAMETERS SETTING

clearvars -except V U BV D F1G
numIterations = 10;
designvars = 'counts';
% 1) 'theta' 2) 'pas' (phase amp midpoint) 3) 'counts'
classes = 'lick';
% 1) 'gonogo' 2) 'FAvsCR' 3) 'lick' 4) allBehavTypes
normalization = 'none';
% 1) 'whiten' 2) 'none';
sample ='random';
% 1) 'bias' (takes 70% from each class for train) 2) 'random' just takes
% random 70% to train
balance = 'on';
% Design Matrix


accprop=cell(1,length(V));

for rec = 1:length(V)
    switch designvars
        case 'theta'
            hx = [V(rec).var.hit{1}'];
            FAx = [V(rec).var.FA{1}'];
            CRx = [V(rec).var.CR{1}'];
            hy = ones(size(hx,1),1);
            FAy = ones(size(FAx,1),1)+1;
            CRy1 = ones(size(CRx,1),1)+1;
            
            F1color = 'b';
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
            
        case 'counts'
                
                hx = [V(rec).touchNum.hit'];
                hy = ones(size(hx,1),1);
                FAx = [V(rec).touchNum.FA'];
                FAy = ones(size(FAx,1),1)+1;
                CRx = [V(rec).touchNum.CR'];
                CRy1 = ones(size(CRx,1),1)+1;
                F1color = 'r';

        case 'ubered'
            Ttmp=[V(rec).touchNum.hit];
            noTouches = find(Ttmp==0);
            tmp2 = [Ttmp zeros(1,length(Ttmp))];            
            for i = 1:length(Ttmp)
                tmp0 = [zeros(1,i) Ttmp zeros(1,length(Ttmp)-i)];
                tmp2 = [tmp2 ;tmp0];
            end   
            tIdx = sum(tmp2(:,1:length(Ttmp)));
            theta = V(rec).var.hit{1}
            lastTheta = theta(tIdx);
            
    end
    
    %     [FAx,FAy,CRx,CRy1] = FACRBalance([FAx ;hx],CRx);
    
    switch classes
        case 'gonogo'
            DmatX=[hx;FAx;CRx]; DmatY=[hy;FAy;CRy1]; %complete design matrix for govsnogo
        colors = {'b','r'};
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
                [FAx,FAy,CRx,CRy1] = FACRBalance([hx ; FAx],CRx);
                DmatX = [FAx;CRx]; DmatY = [FAy-1;CRy1];
            else
                DmatX = [hx;FAx;CRx]; DmatY = [hy;FAy-1;CRy1];
            end
            colors = {'b','r'};
            g1 = [hy;FAy-1];g2 = [CRy1];
        case 'allBehavTypes'
            DmatX=[hx;FAx;CRx]; DmatY = [hy;FAy;CRy1+1];
    end
    
    switch normalization
        case 'whiten'
            DmatX = filex_whiten(DmatX);
    end
    
    
    
    clear Acc
    clear opt_thresh
    
    for reit = 1:numIterations
        rando = randperm(length(DmatX));
        tmpDmatX=DmatX(rando,:);tmpDmatY=DmatY(rando,:);
        
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
                
                [pred,opt_thresh(reit)]=ML_predictOneVsAll(thetas,tmpDmatX(test,:)...
                    ,tmpDmatY(test,:),'Max');
                Acc(reit)= mean(double(pred == tmpDmatY(test))) * 100;
                F1s(reit,:) = F1score(pred,tmpDmatY(test),2);
                
                
                accprop{rec}=[accprop{rec} ; pred tmpDmatY(test)];
                
            case 'random'
                % use this for GO vs NOGO since it doesnt bias sampling
                [thetas,cost,~] = ML_oneVsAll(tmpDmatX(1:end*.7,:),tmpDmatY(1:end*.7,:),numel(unique(DmatY)),0);
                Bfeat{rec}.theta{reit}=thetas;
                [pred,opt_thresh(reit)]=ML_predictOneVsAll(thetas,tmpDmatX(end*.7:end,:)...
                    ,tmpDmatY(end*.7:end,:),'Max');
                Acc(reit)= mean(double(pred == tmpDmatY(end*.7:end) )) * 100;
                F1s(reit,:) = F1score(pred,tmpDmatY(end*.7:end),2);
                
                accprop{rec}=[accprop{rec} ; pred tmpDmatY(end*.7:end)];
                
        end
        %fprintf('\nTraining Set Accuracy: %f\n', mean(double(pred == tmpDmatY(end*.6:end))) * 100)
        
    end
    
    train_F1s(rec,:) = nanmean(nansum(F1s,2));
    trainF1sstd(rec,:)=nanstd(nansum(F1s,2));
    train_predOpt(rec)=mean(opt_thresh);
    train_Acc(rec) = mean(Acc);
    train_std(rec) = std(Acc);
    
    
end


% if strcmp(designvars,'theta')
%     F1G{1} = train_F1s ;
% elseif strcmp(designvars,'counts')
%     F1G{2} = train_F1s ;
%     figure(20);hold on;scatter(F1G{1},F1G{2},'ro');
%     set(gca,'xlim',[0 2],'ylim',[0 2],'ytick',0:.5:2);
%     xlabel('F1 Score Theta');ylabel('F1 Score Touch Count');
%     hold on; plot([0 2],[0 2],'-.k')
% %     legend('Continuous','Discrete','location','northwest')
%     
% end

figure(22);hold on;
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
hold on; scatter(1:length(V),predprop(:,1),colors{1});
scatter(1:length(V),predprop(:,2),colors{2})
legend('Model Accuracy','Go Prediction Accuracy','No Go Prediction Accuracy','location','southeast')
title([U{rec}.meta.layer ' ' designvars ' ' classes])
print(figure(22),'-dtiff',['Z:\Users\Jon\Projects\Characterization\' U{rec}.meta.layer '\Figures\'  U{rec}.meta.layer '_' classes '_' designvars])

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
            
            
            figure(12);subplot(3,3,rec)
            h2a=plot(x,firstvar/sum(firstvar),colors{1});h2a.Color(4)=0.5;
            hold on;h2a=plot(x,secondvar/sum(secondvar),colors{2});h2a.Color(4)=0.5;
            
            for db = [1: 2]
                ms=cell2mat(Bfeat{rec}.theta);
                coords=mean(reshape(ms(db,:)',2,numIterations),2);
                y= (exp(coords(1)+coords(2)*x)) ./ (1 + exp(coords(1)+coords(2)*x))  ;
                y1=train_predOpt(rec);
                hold on; plot(x,y,['-.o' colors{db}]);
           
                set(gca,'xlim',[min(x) max(x)],'ylim',[0 1]);
                xlabel('Counts');
            end
            if rec ==3
                legend('Lick', 'No Lick','Lick Classifier','No Lick Classifier')
            elseif rec == 4
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
            
            figure(10);%hold on;subplot(2,3,rec);
            scatter3(alldat(1:length(hx),1),alldat(1:length(hx),2),alldat(1:length(hx),3),'b')
            hold on;scatter3(alldat(length(hx):end,1),alldat(length(hx):end,2),alldat(length(hx):end,3),'r')
            %             hold on;scatter3(centroid(1,1),centroid(1,2),centroid(1,3),'k','linewidth',10)
            %             hold on;scatter3(centroid(2,1),centroid(2,2),centroid(2,3),'b','linewidth',10)
            
            plot_x = [min(alldat(:,1))-2, min(alldat(:,1))-2, max(alldat(:,1))+2, max(alldat(:,1))+2]; %ranges for amplitude
            plot_z = [-3 ,3,3,-3];
            plot_y = (-1/coords(3)) .* (coords(1) + (coords(2).*plot_x) + (coords(4).*plot_z) - log(train_predOpt(rec)/(1-train_predOpt(rec)))); % log p(go trial) to calculate decision boundary
            
            hold on; fill3(plot_x, plot_y, plot_z,'k'); alpha (.4)
            xlabel('Amplitude');ylabel('Midpoint');zlabel('Phase')
            
            
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
%% Plotting most predictive motor feature
figure(92);filex_barwitherr(optfeatstd',optfeat')
xlabel('Mouse Number');title('Most Predictive Motor Feature');
ylabel('Weight');
legend('Bias','Amplitude','Midpoint','Phase','Location','northeast')
%% FA:CR DISTRIBUTION
% myC=[0 1 0;1 0 0]; % make a colors list
% H=bar(FCratio, 'stack');
% for k=1:2
%     set(H(k),'facecolor',myC(k,:))
% end
% xlabel('Mouse Number'); ylabel('Proportion of Trials')
