%change varname/var to variable you want to plot all of
%change plotrow/plotcol


varname='thetaall';
var=thetaall;

figure(39);clf;
plotrow=8;
plotcol=8;
%%
if strcmp(varname,'ampall')|strcmp(varname,'spall')|strcmp(varname,'thetaall')
    for p=1:length(U)
        subplot(plotrow,plotcol,p);
        imagesc(var{p}(:,2:77))
        hold on
        colormap(gca,parula)
        caxis([0 max(var{p}(:,79))])
        if strcmp(varname,'thetaall')
            range = var{p}(:,1);
        else
            range = round(var{p}(:,1));
        end
        plot([25 25],[25 0],'w:')
        set(gca,'Ydir','normal','ytick',(1:length(var{p}(:,1))),'yticklabel',[range],'xtick',(0:25:75),'xticklabel',(-25:25:50));
        for k=1:size(var{p},1)
            text(20,k,num2str(var{p}(k,78)),'FontSize',8,'Color','white')
        end
        axis('square')
    end
else
    for p=1:length(U)
        subplot(plotrow,plotcol,p);
        imagesc(phaseall{p}(:,2:77))
        hold on
        plot([25 25],[25 0],'w:')
        set(gca,'xtick',[0:5:75]);
        set(gca,'xticklabel',[-25:5:50]);
        title('Phase at Touch')
        colormap(gca,parula)
        caxis([0 max(phaseall{p}(:,79))])
        set(gca,'Ydir','normal','ytick',[1:2.5:12],'yticklabel',{'-pi','-pi/2',0,'pi/2','pi'},'xtick',(0:25:75),'xticklabel',(-25:25:50))
        for k=1:size(phaseall{p},1)
            text(20,k+.1,num2str(phaseall{p}(k,78)),'FontSize',8,'Color','white')
        end
        axis('square')
    end
end