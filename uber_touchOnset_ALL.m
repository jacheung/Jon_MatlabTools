figure(39);clf;
plotrow=5;
plotcol=7;   

for k = 1:length(touchONN)
    subplot(plotrow,plotcol,k);
    bar(-50:150,sum(touchONN{k})/size(touchONN{k},1),'k');
    set(gca,'xlim',[-50 150]);
    if k==1
    xlabel('Time from all touch onset (ms)')
    ylabel('spks / ms')
    end
end    
