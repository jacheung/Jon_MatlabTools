%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rebuild psychometric curves using U array (real) and for model. Need to
% input U array and PredArray(1st column = motorPos, and 2nd column =
% prediction)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [RMSEgroup] = rebuildPsychoflip(U,V,predArray)



for rec = 1:length(U)
    
    if strcmp(U{rec}.meta.layer,'D')
        nolickmean = predArray{rec}(predArray{rec}(:,1) == -1,2);
        lickmean = predArray{rec}(predArray{rec}(:,1) == 1,2);
        gos = U{rec}.meta.trialType ==1;
        nogos = U{rec}.meta.trialType == 0 ;
        
        modelmean = [mean(nolickmean) mean(lickmean)];
        modelstd = [std(nolickmean) std(lickmean)];
        
        realmean = [1-mean(U{rec}.meta.trialCorrect(nogos)) mean(U{rec}.meta.trialCorrect(gos))];
        
        RMSE = sqrt(nanmean((realmean-modelmean).^2));
        RMSEdec = sqrt((realmean-modelmean).^2);
        
        figure(5);subplot(2,5,rec)
        scatter([-1 1],realmean,'filled','r')
        hold on;errorbar([-1 1],modelmean,modelstd,'ko')
        set(gca,'xlim',[-1.5 1.5],'ylim',[0 1])
        
        if rec ==5
            legend('Mouse Performance','Model Performance','location','southeast')
        end
        
        text(.5,.4,['RMSE = ' num2str(RMSE)])
        
    elseif strcmp(U{1}.meta.layer,'SM')
        
        [sorted]= binslin(predArray{rec}(:,1),predArray{rec},'equalE',12,U{rec}.meta.ranges(1),U{rec}.meta.ranges(2));
        lickmean=cell2mat(cellfun(@mean,sorted,'uniformoutput',0));
        lickstd=cell2mat(cellfun(@std,sorted,'uniformoutput',0));
        xranges = U{rec}.meta.ranges(1):10000:U{rec}.meta.ranges(2);
        
        real=[U{rec}.meta.motorPosition;V(rec).trialNums.matrix(5,:)]';%only taking lick row and motorPos row
        [realsorted]= binslin(real(:,1),real,'equalE',12,U{rec}.meta.ranges(1),U{rec}.meta.ranges(2));
        reallickmean=cell2mat(cellfun(@(x) mean(x,1),realsorted,'uniformoutput',0));
        reallickstd=cell2mat(cellfun(@(x) std(x,0,1),realsorted,'uniformoutput',0));
        
        RMSE = sqrt(nanmean((reallickmean(:,2)-lickmean(:,2)).^2));
        RMSEdec = sqrt((reallickmean(:,2)-lickmean(:,2)).^2);
        
        figure(5);subplot(2,5,rec)
        xranges(isnan(lickmean(:,1)))=[];
        lickmean(isnan(lickmean(:,1)),:)=[];
        lickstd(isnan(lickstd(:,1)),:)=[];
        reallickmean(isnan(reallickmean(:,1)),:)=[];
        scatter(xranges(1),flipud(lickmean(1,2)),'ko','filled')
        hold on;filex_shadedErrorBar(xranges(end-(length(lickmean)-2):end),flipud(lickmean(2:end,2)),flipud(lickstd(2:end,2)),'k')
        scatter(xranges(1),flipud(reallickmean(1,2)),'ro','linewidth',1)
        plot(xranges(end-(length(lickmean)-2):end),flipud(reallickmean(2:end,2)),'r','linewidth',2)
        
        if rec ==5
            legend('Mouse Performance','Model Performance','location','southeast')
        end
        
        
        set(gca,'xtick',[xranges(1) xranges(end-(length(lickmean)-2)) xranges(end)],'xticklabel',[-1 0 1],'ylim',[0 1],'xlim',[xranges(1)-5000 xranges(end)+5000])
        xlabel('Motor Position')
        ylabel('Lick Probability')
        text(xranges(end)*.5,.4,['RMSE = ' num2str(RMSE)])
        
    elseif strcmp(U{1}.meta.layer,'BV')
        [sorted]= binslin(predArray{rec}(:,1),predArray{rec},'equalE',12,U{rec}.meta.ranges(1),U{rec}.meta.ranges(2));
        lickmean=cell2mat(cellfun(@mean,sorted,'uniformoutput',0));
        lickstd=cell2mat(cellfun(@std,sorted,'uniformoutput',0));
        xranges = U{rec}.meta.ranges(1):10000:U{rec}.meta.ranges(2);
        
        real=[U{rec}.meta.motorPosition;V(rec).trialNums.matrix(5,:)]';%only taking lick row and motorPos row
        [realsorted]= binslin(real(:,1),real,'equalE',12,U{rec}.meta.ranges(1),U{rec}.meta.ranges(2));
        reallickmean=cell2mat(cellfun(@(x) mean(x,1),realsorted,'uniformoutput',0));
        reallickstd=cell2mat(cellfun(@(x) std(x,0,1),realsorted,'uniformoutput',0));
        RMSE = sqrt(nanmean((reallickmean(:,2)-lickmean(:,2)).^2));
        RMSEdec = sqrt((reallickmean(:,2)-lickmean(:,2)).^2);
        
        
        figure(5);subplot(2,5,rec)
        plot(xranges',flipud(reallickmean(:,2)),'r','linewidth',5)
        hold on;filex_shadedErrorBar(xranges, flipud(lickmean(:,2)),flipud(lickstd(:,2)),'k');
        hold on; plot(xranges',flipud(reallickmean(:,2)),'r','linewidth',5)
        set(gca,'xtick',[xranges(1) xranges(6) xranges(11)],'xticklabel',[-1 0 1],'ylim',[0 1],'xlim',[xranges(1)-5000 xranges(end)+5000],'ytick',[0 .5 1])
        if rec ==5
            legend('Mouse Performance','Model Performance','location','southeast')
%         elseif rec == 3
%             title([U{rec}.meta.layer ' ' designvars ' ' classes])
%             
        end
        
%         xlabel('Motor Position')
%         ylabel('Lick Probability')
        text(xranges(end)*.5,.3,['RMSE = ' num2str(RMSE)])
        
    end
    
    RMSEgroup(rec) = RMSE;
    
%     if strcmp(designvars,'counts')
%         RMSEgroup(rec,1) = RMSE;
% %         RMSEdecomp{1}(:,rec) = RMSEdec;
%     elseif strcmp(designvars,'theta')
%         RMSEgroup(rec,2) = RMSE;
% %         RMSEdecomp{2}(:,rec) = RMSEdec;
%     elseif strcmp(designvars,'ubered')
%         RMSEgroup(rec,3) = RMSE;
% %         RMSEdecomp{3}(:,rec) = RMSEdec;
%     end
%     
%     if size(RMSEgroup,2)==3
% 
%         if strcmp(U{1}.meta.layer,'BV')
%             R{3} = RMSEgroup;
%         elseif strcmp(U{1}.meta.layer,'SM')
%             R{2} = RMSEgroup;
%         elseif strcmp(U{1}.meta.layer,'D')
%             R{1} = RMSEgroup;
%         end
%     end
    
    
end
set(gcf, 'Units', 'pixels', 'Position', [0, 0, 2000, 1000]);
