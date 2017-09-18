% load('Z:\Users\Jon\Projects\Characterization\BV')
% U=BV;
%
countsALL=cell(1,length(U));
gocorrALL=cell(1,length(U));
nogocorrALL=cell(1,length(U));
godurALL=cell(1,length(U));
nogodurALL=cell(1,length(U));

for rec=2
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
         print(figure(70+f),'-dtiff',['Z:\Users\Jon\Projects\Characterization\' U{rec}.meta.layer '\Figures\' U{rec}.meta.layer num2str(rec) '_' V(rec).varNames{f}])
  
    end
         
    [counts, gocorr, nogocorr,godur,nogodur] = assist_countsNduration(U{rec},prelixGo,prelixNoGo,duration);
    print(figure(123),'-dtiff',['Z:\Users\Jon\Projects\Characterization\' U{rec}.meta.layer '\Figures\' U{rec}.meta.layer num2str(rec) '_countsNduration' ])
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





