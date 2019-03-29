
type = {SM,BV}
for d = 2
    U = type{d};
    
    
    for k = 1:length(U)
        array = U{k};
        [tperWhisksorted, trialNums,stats] = uber_whisksThatTouch(array);
        accountedTouches{d}(k) = stats.proportionTouches;
        colors = {'b','g','r','k'};
        
        goprops = [ tperWhisksorted{1} ;tperWhisksorted{2}];
        nogoprops = [ tperWhisksorted{4} ;tperWhisksorted{3}];
        
            gotmp = histc(goprops,0:.1:1);
            nogotmp = histc(nogoprops,0:.1:1);
            
            figure(480);subplot(2,5,k)
            bar(0:.1:1,gotmp,'g');alpha(.5)
            hold on;
            bar(0:.1:1,nogotmp,'k');alpha(.5)
        
            set(gca,'xlim',[-.05 1],'xtick',0:.5:1)
        
        
    end
    
    figure(830);subplot(2,1,d)
    bar(1:length(U),accountedTouches{d},'k')
    set(gca,'ylim',[0 1],'ytick',[0:.5:1])
    pause
end

