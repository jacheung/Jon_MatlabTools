%% tuning quant using cdfs and firing rate (need to consider choosing optimal window to count spikes).

CDFquant=zeros(length(pop),1);
for d = 1:length(pop)
binnedrates=mean(pop{d}.theta.spikes(:,51:76),2)*1000;


    current=sort(binnedrates./sum(binnedrates));
    
    previousVal = 0;
    cdfvals=zeros(length(current),1);
    for k = 1:length(current)
        previousVal = previousVal+current(k);
        cdfvals(k)=previousVal;
    end  

cdfvals = [0; cdfvals];

    figure(4030);clf;
    plot(0:length(cdfvals)-1,cdfvals)
    set(gca,'xlim',[0 length(cdfvals)-1],'ylim',[0 1])
    hold on; plot([0 length(cdfvals)-1],[0 1],'-.k')
    title(['n=' num2str(d) ' ' num2str(sum(cdfvals))])
    
    CDFquant(d)= sum(cdfvals);
end 

 [~,idx] = sort(CDFquant);
 

figure(32);clf
fields = [    {'theta'}  ];
for d= 1:length(pop)
    currCell=idx(d);
    
    
        figure(32);subplot(5,8,d);
       imagesc(imgaussfilt(pop{currCell}.(fields{1}).spikes,gaussFilt,'padding','replicate'));
    colormap(gca,parula);
    set(gca,'Ydir','normal','ytick',(1:length(pop{currCell}.(fields{1}).range)),'yticklabel',[pop{currCell}.(fields{1}).range],...
        'xtick',(0:25:length(window)),'xticklabel',[min(window):25:max(window)],'xlim',[0 length(window)]);
    for k=1:size(pop{d}.(fields{1}).counts,1)
        text(20,k,num2str(pop{d}.(fields{1}).counts(k)),'FontSize',8,'Color','white')
    end
    hold on;plot([sum(window<0) sum(window<0)],[length(pop{d}.(fields{1}).range) 0],'w:')

end