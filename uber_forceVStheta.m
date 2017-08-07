load('Z:\Users\Jon\Projects\Characterization\BV')
U=BV;
%
countsALL=cell(1,length(U));
gocorrALL=cell(1,length(U));
nogocorrALL=cell(1,length(U));
godurALL=cell(1,length(U));
nogodurALL=cell(1,length(U));

for rec=1:length(U)
    varcell=cell(1,3);
    prevarcell=varcell;
    %[varspikes] = assist_varAtTouch(U{rec},-25:50);
    var=[1 2 6];
    for f=1:length(var)
        [varForce,prevarForce] = assist_varAtTouchFORCE(U{rec},-25:50,var(f));
        %varcell{var}=varspikes; varcell{2}=varForce;
        varcell{f}=varForce;
        prevarcell{f}=prevarForce;
        
    end
    varcell=[varcell prevarcell];
    figure(39); clf
    normbinmin=10; %min in each bin to use for normalizing
    normbinwin=30; %window for normalizing from 0ms:binwin
    
    xtmp=cell(1,length(varcell));
    
    for d=1:numel(varcell)
        %[sorted sortedBy binBounds]=binslin(varcell{d}(:,1),varcell{d}(:,7:82),'equalE',41,-100,100);
        [sorted sortedBy binBounds]=binslin(varcell{d}(:,1),varcell{d}(:,2:end),'equalE',41,-100,100);
        
        thetarange=[-97.5:5:97.5];
        
        thetaresponse=zeros(length(sorted),79);
        for j=1:length(sorted)
            thetaresponse(j,2:77)=mean(sorted{j},1);
            thetaresponse(j,78)=size(sorted{j},1);
        end
        thetaresponse(:,1)=thetarange;
        %thetatouchspks = thetatouchspks(~any(isnan(thetatouchspks),2),:); %remove NaN rows
        thetaresponse = thetaresponse(~any(isnan(thetaresponse(:,2)),2),:); %remove NaN rows USING column for theta
        
        for k=1:size(thetaresponse,1)
            thetaresponse(k,2:77)=smooth(thetaresponse(k,2:77));
            if thetaresponse(k,78)>normbinmin
                thetaresponse(k,79)=max(thetaresponse(k,26:26+normbinwin)); %normalizing window
            end
        end
        
        figure(39);h1 = subplot(2,3,d);
        imagesc(thetaresponse(:,2:77))
        hold on
        colormap(gca,parula)
        hCBar=colorbar;
        caxis([0 max(thetaresponse(:,79))])
        plot([25 25],[25 0],'w:')
        set(gca,'Ydir','normal','ytick',(1:length(thetaresponse(:,1))),'yticklabel',[thetaresponse(:,1)],'xtick',(0:25:75),'xticklabel',(-25:25:50));
        for k=1:size(thetaresponse,1)
            text(20,k,num2str(thetaresponse(k,78)),'FontSize',8,'Color','white')
        end
        axis('square')
        if d==2
            title('DTheta       pretouch vel      DKappa')
            xlabel('time from touch onset')
        end
    end
    V(rec).varNames = {'theta','velocity','amplitude','setpoint','phase','deltaKappa'};
    array=U{rec};
    [~ ,prelixGo, prelixNoGo, ~ ,duration ,~] = assist_predecisionVar(array);
    varx=[1:6];
    for f=1:length(varx)
        [~, ~,V(rec).hit{f}, V(rec).miss{f}, V(rec).FA{f}, V(rec).CR{f},~,~] = assist_vardistribution(U{rec},varx(f),prelixGo,prelixNoGo,[-25:0],[5:25],'on');
        if f == 1
        V(rec).touchNum.hit=cellfun(@numel,V(rec).hit{1}); % count number of predecision touches using theta vals
        V(rec).touchNum.miss=cellfun(@numel,V(rec).miss{1});
        V(rec).touchNum.FA=cellfun(@numel,V(rec).FA{1});
        V(rec).touchNum.CR=cellfun(@numel,V(rec).CR{1});
        end
        V(rec).hit{f}=cell2mat(V(rec).hit{f}); %predecision variables for hit trials
        V(rec).miss{f} = cell2mat(V(rec).miss{f});
        V(rec).FA{f} = cell2mat(V(rec).FA{f});
        V(rec).CR{f} = cell2mat(V(rec).CR{f});
        
                 print(figure(70+f),'-dtiff',['Z:\Users\Jon\Projects\Characterization\BV\Figures\BV' num2str(rec) '_' U{rec}.varNames{f}])
    end
    [counts, gocorr, nogocorr,godur,nogodur] = assist_countsNduration(U{rec},prelixGo,prelixNoGo,duration);
    %print(figure(123),'-dtiff',['Z:\Users\Jon\Projects\Characterization\' U{rec}.meta.layer{1} '\Figures\' U{rec}.meta.layer{1} num2str(rec) '_countsNduration' ])
    countsALL{rec}=counts;
    gocorrALL{rec}=gocorr;nogocorrALL{rec}=nogocorr; % 2 columns, 1 = count num, 2 = p(lick)
    godurALL{rec}=godur;nogodurALL{rec}=nogodur; %3 rows, 1 = ?, 2 = duration, 3 = p(lick)
    
    
    hx = [V(rec).touchNum.hit'];
    FAx = [V(rec).touchNum.FA'];
    CRx = [V(rec).touchNum.CR'];
    gotouch = histc(hx,0:10);nogotouch = histc([FAx;CRx],0:10);
    figure(40);clf;plot(0:10,gotouch/sum(gotouch),'b');hold on;plot(0:10,nogotouch/sum(nogotouch),'r')
    xlabel('Number predecision touches');ylabel('Proportion of trials');title('Predecision Touches and Go/Nogo')
    
end

%% Design Matrix
clearvars -except V

for rec = 1:length(V)
    
     hx = [V(rec).hit{1}']; 
%     ntmp=find(V(rec).hit{5}<=0);ptmp=find(V(rec).hit{5}>0);
%     hx = [V(rec).hit{3}(ntmp)' V(rec).hit{4}(ntmp)' V(rec).hit{5}(ntmp)'];
    hy = ones(size(hx,1),1);
    
    FAx = [V(rec).FA{1}'];
%     Fntmp=find(V(rec).FA{5}<=0);Fptmp=find(V(rec).FA{5}>0);
%     FAx = [V(rec).FA{3}(Fntmp)' V(rec).FA{4}(Fntmp)' V(rec).FA{5}(Fntmp)'];
    FAy = ones(size(FAx,1),1)+1;

    CRx = [V(rec).CR{1}']; 
%     Cntmp=find(V(rec).CR{5}<=0);Cptmp=find(V(rec).CR{5}>0);
%     CRx = [V(rec).CR{3}(Cntmp)' V(rec).CR{4}(Cntmp)' V(rec).CR{5}(Cntmp)']; 
    CRy1 = ones(size(CRx,1),1)+1;
    %FOR FA VS CR
     %RESAMPLING so CR and FA will be balanced
    [FAx,FAy,CRx,CRy1] = FACRBalance(FAx,CRx);
    
    %FOR COUNTS 
%     hx = [V(rec).touchNum.hit'];
%     hy = ones(size(hx,1),1); 
%     FAx = [V(rec).touchNum.FA'];
%     FAy = ones(size(FAx,1),1)+1;
%     CRx = [V(rec).touchNum.CR'];
%     CRy1 = ones(size(CRx,1),1)+1;
%     
%     [FAx,FAy,CRx,CRy1] = FACRBalance(FAx,CRx);
    


    
    %COMPLETE DESIGN MATRIX 
%     DmatX{1}=filex_whiten([hx;FAx;CRx]); DmatY{1}=[hy;FAy;CRy1]; %complete design matrix for govsnogo
    % DmatX{2}=[hx;FAx;CRx]; DmatY{2}=[hy;FAy;CRy]; %complete design matrix for trial type discrim
    %DmatX{1}=filex_whiten([FAx;CRx]); DmatY{1} = [FAy-1;CRy1]; %complete design matrix for FA CR
    
%      DmatX{1} = ([hx;FAx;CRx]) ;DmatY{1}=[hy;FAy;CRy1]; %complete design matrix for touch counts go/nogo
    DmatX{1}=([FAx;CRx]); DmatY{1} = [FAy-1;CRy1]; %complete design matrix for touch count FA CR
    
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

            
            %fprintf('\nTraining Set Accuracy: %f\n', mean(double(pred == tmpDmatY(end*.6:end))) * 100)
            
        end
        
        train_predOpt{z}(rec)=mean(opt_thresh);
        train_Acc{z}(rec) = mean(Acc);
        train_std{z}(rec) = std(Acc);
    end
    
end

figure(22);
errorbar(1:length(V),train_Acc{1},train_std{1},'bo')
% hold on; errorbar(1:length(V),train_Acc{1},train_std{1},'go')
% hold on; errorbar(1:length(V),train_Acc{3},train_std{3},'ro')
%title('Prediction Accuracy of FA vs CR')
xlabel('Mouse Number');ylabel('% Accuracy')
set(gca,'ylim',[25 100])
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
            hx = [V(rec).hit{1}'];
            FAx = [V(rec).FA{1}'];
            CRx = [V(rec).CR{1}'];
            
            ms=cell2mat(Bfeat{rec}.theta{1});
            coords=mean(reshape(ms(1,:)',2,10),2);
            x = [min([hx;FAx;CRx])-2:1:max([hx;FAx;CRx])+2];
            y= (exp(coords(1)+coords(2)*x)) ./ (1 + exp(coords(1)+coords(2)*x))  ;
            y1=train_predOpt{1}(rec);
            firstvar = histc(FAx,x);secondvar = histc([CRx],x);
            
            figure(12);subplot(2,3,rec);plot(x,firstvar/sum(firstvar),'g')
            hold on;plot(x,secondvar/sum(secondvar),'r')
            hold on; plot(x,y,'k');
            hold on; plot([min(x) max(x)],[y1 y1],'-.ok')
            set(gca,'xlim',[min(x) max(x)],'ylim',[0 1]);
            xlabel('Theta at Touch');ylabel('Proportion of Trials AND P(FA) trial')
            if rec ==3 
            legend('FA','CR','FA Classifier','Decision Boundary')
            end
            
            % FOR TOUCH COUNTS
%             ms=cell2mat(Bfeat{rec}.theta{1});
%              coords=mean(reshape(ms(1,:)',2,10),2);
%             hx = [V(rec).touchNum.hit'];
%             FAx = [V(rec).touchNum.FA'];
%             CRx = [V(rec).touchNum.CR'];
%               x = [0:1:10];
%             firstvar = histc(FAx,x);secondvar = histc([CRx],x);
%             y= (exp(coords(1)+coords(2)*x)) ./ (1 + exp(coords(1)+coords(2)*x))  ;
%             y1=train_predOpt{1}(rec);
%             
%             figure(12);subplot(2,3,rec);plot(0:10,firstvar/sum(firstvar),'g')
%             hold on;plot(0:10,secondvar/sum(secondvar),'r')
%             hold on; plot(x,y,'k');
%             hold on; plot([0 10],[y1 y1],'-.ok')
%             set(gca,'xlim',[0,10],'ylim',[0 1]);
%             xlabel('Number of Touches');ylabel('Proportion of Trials AND P(FA) trial')
%             if rec ==3 
%             legend('FA','CR','FA Classifier','Decision Boundary')
%             end
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

    end


end
%% FA:CR DISTRIBUTION 
myC=[0 1 0;1 0 0]; % make a colors list
H=bar(FCratio, 'stack');
for k=1:2
    set(H(k),'facecolor',myC(k,:))
end
xlabel('Mouse Number'); ylabel('Proportion of Trials')

%% Population Plotting
% Counts Population
xcountsALL=cell2mat(countsALL);
gocounts=mean(xcountsALL(:,1:2:end),2);
nogocounts=mean(xcountsALL(:,2:2:end),2);
allcounts=[gocounts nogocounts];
allcounts=[allcounts(1:10,:);nanmean(allcounts(11:end,:),1)];%condensing anything >10 touches into a single variable

gocorr=cell2mat(gocorrALL);nogocorr=cell2mat(nogocorrALL);
gocorrmean=nanmean(gocorr(:,2:2:end),2);nogocorrmean=nanmean(nogocorr(:,2:2:end),2);
corrmeans=[(0:19)' gocorrmean nogocorrmean];
corrmeans=[corrmeans(1:10,:);nanmean(corrmeans(11:end,:))]; corrmeans(end,1)=10;%condensing >10 pre decision to just one single value

figure(124);clf;subplot(1,2,1)
bar(0:10,[allcounts])
hold on; plot(corrmeans(:,1),corrmeans(:,2),'-o') %plot go pLick
hold on; plot(corrmeans(:,1),1-corrmeans(:,3),'r-o')%plot nogo pLick
xlabel('Number of contacts per trial');ylabel('P(Lick) or Touches per Trials')
title('Population Predecision Contacts x P(lick)')
set(gca,'xlim',[0 10])

% Duration Population
xgodurALL=cell2mat(godurALL');xnogodurALL=cell2mat(nogodurALL');
godurPOP=nanmean(xgodurALL(1:3:end,:));nogodurPOP=nanmean(xnogodurALL(1:3:end,:));
godurcorrPOP=nanmean(xgodurALL(3:3:end,:));nogodurcorrPOP=nanmean(xnogodurALL(3:3:end,:));
godurcorrPOP=[godurcorrPOP(1:8) nanmean(godurcorrPOP(1,9:end))];
nogodurcorrPOP=1-[nogodurcorrPOP(1:8) nanmean(nogodurcorrPOP(1,9:end))];%add -1 since looking at P(lick)

durPOP=[godurPOP;nogodurPOP]';%condensing touches >50ms into a single variable
durPOP=[durPOP(1:8,:) ;nanmean(durPOP(9:end,:))];


figure(124);subplot(1,2,2);
bar(godurALL{1}(2,1:9),durPOP)
hold on; plot(godurALL{1}(2,1:9),godurcorrPOP,'-o')
hold on; plot(godurALL{1}(2,1:9),nogodurcorrPOP,'-ro')
xlabel('Duration of Contacts (ms)');ylabel('P(Lick) or Touches per Trial')
title('Population Duration Contacts x P(lick)')
set(gca,'xtick',[0:10:55],'xlim',[0 35])





