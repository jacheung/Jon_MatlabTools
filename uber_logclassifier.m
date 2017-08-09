clear all; close all;clc;

load('Z:\Users\Jon\Projects\Characterization\BV')
U=BV;
[V] = classifierWrapper(U);

%% Design Matrix
clearvars -except V U

for i = 1:length(U)
    pc(i)=mean(U{i}.meta.trialCorrect);
end
accprop=cell(1,length(V))
for rec = 1:length(V)
    
    %      hx = [V(rec).hit{1}'];
    %     ntmp=find(V(rec).hit{5}<=0);ptmp=find(V(rec).hit{5}>0);
    % %     hx = [V(rec).hit{3}(ntmp)' V(rec).hit{4}(ntmp)' V(rec).hit{5}(ntmp)'];
    %     hy = ones(size(hx,1),1);
    %
    %     FAx = [V(rec).FA{1}'];
    %     Fntmp=find(V(rec).FA{5}<=0);Fptmp=find(V(rec).FA{5}>0);
    % %     FAx = [V(rec).FA{3}(Fntmp)' V(rec).FA{4}(Fntmp)' V(rec).FA{5}(Fntmp)'];
    %     FAy = ones(size(FAx,1),1)+1;
    %
    %     CRx = [V(rec).CR{1}'];
    %     Cntmp=find(V(rec).CR{5}<=0);Cptmp=find(V(rec).CR{5}>0);
    % %     CRx = [V(rec).CR{3}(Cntmp)' V(rec).CR{4}(Cntmp)' V(rec).CR{5}(Cntmp)'];
    %     CRy1 = ones(size(CRx,1),1)+1;
    %     %FOR FA VS CR
    %      %RESAMPLING so CR and FA will be balanced
    %     [FAx,FAy,CRx,CRy1] = FACRBalance(FAx,CRx);
    %
    %     %FOR COUNTS
    hx = [V(rec).touchNum.hit'];
    hy = ones(size(hx,1),1);
    FAx = [V(rec).touchNum.FA'];
    FAy = ones(size(FAx,1),1)+1;
    CRx = [V(rec).touchNum.CR'];
    CRy1 = ones(size(CRx,1),1)+1;
    
    [FAx,FAy,CRx,CRy1] = FACRBalance(FAx,CRx);
    %
    
    
    
    %COMPLETE DESIGN MATRIX
    %     DmatX{1}=filex_whiten([hx;FAx;CRx]); DmatY{1}=[hy;FAy;CRy1]; %complete design matrix for govsnogo
    %     DmatX{1}=filex_whiten([hx;FAx;CRx]); DmatY{1}=[hy;FAy;CRy1+1]; %complete design matrix for trial type discrim
    %     DmatX{1}=filex_whiten([FAx;CRx]); DmatY{1} = [FAy-1;CRy1]; %complete design matrix for FA CR
    
    DmatX{1} = ([hx;FAx;CRx]) ;DmatY{1}=[hy;FAy;CRy1]; %complete design matrix for touch counts go/nogo
    %     DmatX{1}=([hx;FAx;CRx]); DmatY{1}=[hy;FAy;CRy1+1]; %complete design matrix for trial type discrim
    
    %     DmatX{1}=([FAx;CRx]); DmatY{1} = [FAy-1;CRy1]; %complete design matrix for touch count FA CR
    
    for z=1:length(DmatX)
        
        clear Acc
        clear DB
        clear opt_thresh
        
        for reit = 1:10
            rando = randperm(length(DmatX{z}));
            tmpDmatX=DmatX{z}(rando,:);tmpDmatY=DmatY{z}(rando,:);
            
            %             %FOR FA VS CR
            %             sample evenly from FA and CR for training set
            %             CReach = round(numel(CRy1)*.7);
            %             FAeach = round(numel(FAy)*.7);
            %             CRs = find(tmpDmatY==2);
            %             FAs = find(tmpDmatY==1);
            %             train=[FAs(1:FAeach);CRs(1:CReach)];
            %             test = [1:length(tmpDmatY)]';
            %             test(train)=[];
            %
            %             [thetas,cost,~] = ML_oneVsAll(tmpDmatX(train,:),tmpDmatY(train,:),numel(unique(DmatY{z})),0);
            %             Bfeat{rec}.theta{z}{reit}=thetas;
            %
            %             [pred,opt_thresh(reit)]=ML_predictOneVsAll(thetas,tmpDmatX(test,:)...
            %                 ,tmpDmatY(test,:),'F1');
            %              Acc(reit)= mean(double(pred == tmpDmatY(test))) * 100
            
            %             % use this for GO vs NOGO since it doesnt bias sampling
            [thetas,cost,~] = ML_oneVsAll(tmpDmatX(1:end*.7,:),tmpDmatY(1:end*.7,:),numel(unique(DmatY{z})),0);
            Bfeat{rec}.theta{z}{reit}=thetas;
            [pred,opt_thresh(reit)]=ML_predictOneVsAll(thetas,tmpDmatX(end*.7:end,:)...
                ,tmpDmatY(end*.7:end,:),'Max');
            Acc(reit)= mean(double(pred == tmpDmatY(end*.7:end) )) * 100
            accprop{rec}=[accprop{rec} ; pred tmpDmatY(end*.7:end)];
            
            %fprintf('\nTraining Set Accuracy: %f\n', mean(double(pred == tmpDmatY(end*.6:end))) * 100)
            
        end
        
        train_predOpt{z}(rec)=mean(opt_thresh);
        train_Acc{z}(rec) = mean(Acc);
        train_std{z}(rec) = std(Acc);
    end
    
end



figure(22);clf;
hold on; errorbar(1:length(V),train_Acc{1},train_std{1},'ko')
% hold on; errorbar(1:length(V),train_Acc{1},train_std{1},'go')
% hold on; errorbar(1:length(V),train_Acc{3},train_std{3},'ro')
%title('Prediction Accuracy of FA vs CR')
xlabel('Mouse Number');ylabel('% Accuracy')
set(gca,'ylim',[0 100])

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
scatter(1:length(V),predprop(:,2),'g')
scatter(1:length(V),predprop(:,3),'r')




%%
params = 1;
FCratio = [];
for rec = 1:6
    ntmp=find(V(rec).hit{5}<=0);ptmp=find(V(rec).hit{5}>0);
    hx = [V(rec).hit{3}(ntmp)' V(rec).hit{4}(ntmp)' V(rec).hit{5}(ntmp)'];
    hy = ones(size(hx,1),1);
    Fntmp=find(V(rec).FA{5}<=0);Fptmp=find(V(rec).FA{5}>0);
    FAx = [V(rec).FA{3}(Fntmp)' V(rec).FA{4}(Fntmp)' V(rec).FA{5}(Fntmp)'];
    FAy = ones(size(FAx,1),1)+1;
    Cntmp=find(V(rec).CR{5}<=0);Cptmp=find(V(rec).CR{5}>0);
    CRx = [V(rec).CR{3}(Cntmp)' V(rec).CR{4}(Cntmp)' V(rec).CR{5}(Cntmp)'];
    CRy1 = ones(size(CRx,1),1)+1;
    
    %FOR FA VS CR
    [FAx,FAy,CRx,CRy1] = FACRBalance(FAx,CRx);
    FCratio(rec,1) = size([FAx],1)/size([FAx; CRx],1);
    FCratio(rec,2) = 1-FCratio(rec,1);
    
    
    switch params
        case 1
            %FOR TOUCH COUNTS
            hx = [V(rec).touchNum.hit'];
            FAx = [V(rec).touchNum.FA'];
            CRx = [V(rec).touchNum.CR'];
            x = [0:1:10];
            %
            %             hx = [V(rec).hit{1}'];
            %             FAx = [V(rec).FA{1}'];
            %             CRx = [V(rec).CR{1}'];
            %             x = [min([hx;FAx;CRx])-2:1:max([hx;FAx;CRx])+2];
            %
            firstvar = histc(hx,x);secondvar = histc([FAx],x);
            %thirdvar = histc([CRx],x);
            
            figure(12);subplot(2,3,rec)
            h2a=plot(x,firstvar/sum(firstvar),'b');h2a.Color(4)=0.5;
            hold on;h2a=plot(x,secondvar/sum(secondvar),'r');h2a.Color(4)=0.5;
            %hold on;h2a=plot(x,thirdvar/sum(thirdvar),'r');h2a.Color(4)=0.5;
            
            for db = [1:2]
                ms=cell2mat(Bfeat{rec}.theta{1});
                coords=mean(reshape(ms(db,:)',2,10),2);
                y= (exp(coords(1)+coords(2)*x)) ./ (1 + exp(coords(1)+coords(2)*x))  ;
                y1=train_predOpt{1}(rec);
                colors = ['b','r'];
                hold on; plot(x,y,['-.o' colors(db)]);
                %             hold on; plot([0 10],[y1 y1],'-.ok')
                set(gca,'xlim',[min(x) max(x)],'ylim',[0 1]);
                xlabel('Touch Count');ylabel('Proportion of Trials AND P(trialType)')
            end
            if rec ==3
                legend('Hit', 'FA','Hit Classifier','FA Classifier', 'CR Classifier')
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
            
            ms=cell2mat(Bfeat{rec}.theta{1});
            coords=mean(reshape(ms(1,:)',4,10),2);
            alldat=filex_whiten([FAx;CRx]);
            centroid=[mean(alldat(1:length(FAx),:));mean(alldat(length(FAx):end,:))];
            
            figure(10);%hold on;subplot(2,3,rec);
            scatter3(alldat(1:length(FAx),1),alldat(1:length(FAx),2),alldat(1:length(FAx),3),'G')
            hold on;scatter3(alldat(length(FAx):end,1),alldat(length(FAx):end,2),alldat(length(FAx):end,3),'r')
            hold on;scatter3(centroid(1,1),centroid(1,2),centroid(1,3),'k','linewidth',10)
            hold on;scatter3(centroid(2,1),centroid(2,2),centroid(2,3),'b','linewidth',10)
            
            plot_x = [min(alldat(:,1))-2, min(alldat(:,1))-2, max(alldat(:,1))+2, max(alldat(:,1))+2]; %ranges for amplitude
            plot_z = [-3 ,3,3,-3];
            plot_y = (-1/coords(3)) .* (coords(1) + (coords(2).*plot_x) + coords(4).*plot_z - log(train_predOpt{1}(rec)/(1-train_predOpt{1}(rec)))); % log p(go trial) to calculate decision boundary
            hold on; fill3(plot_x, plot_y, plot_z,'k'); alpha (.4)
            xlabel('Amplitude');ylabel('Setpoint');zlabel('Phase')
            
            coords=[];
            ms=cell2mat(Bfeat{rec}.theta{1});
            coords(:,rec)=mean(reshape(ms(1,:)',4,10),2);
            
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
            for db = [1]
                ms=cell2mat(Bfeat{rec}.theta{1});
                coords=mean(reshape(ms(db,:)',4,10),2);
                plot_x = [min(alldat(:,1))-2, min(alldat(:,1))-2, max(alldat(:,1))+2, max(alldat(:,1))+2]; %ranges for amplitude
                plot_z = [-3 ,3,3,-3];
                plot_y = (-1/coords(3)) .* (coords(1) + (coords(2).*plot_x) + coords(4).*plot_z - log(train_predOpt{1}(rec)/(1-train_predOpt{1}(rec)))); % log p(go trial) to calculate decision boundary
                hold on; fill3(plot_x, plot_y, plot_z,'k'); alpha (.4)
                xlabel('Amplitude');ylabel('Setpoint');zlabel('Phase')
            end
            
            ms=cell2mat(Bfeat{rec}.theta{1});
            optfeat(:,rec)=mean(reshape(ms(1,:)',4,10),2);
            
            
    end
    
    
end
%% FA:CR DISTRIBUTION
myC=[0 1 0;1 0 0]; % make a colors list
H=bar(FCratio, 'stack');
for k=1:2
    set(H(k),'facecolor',myC(k,:))
end
xlabel('Mouse Number'); ylabel('Proportion of Trials')
