Uraw=U;
[touchCells,peaks] = touchCell(U,2,.5);
selectedCells = find(touchCells==1);
U = U(selectedCells);
%%
figure(329);clf
swin = -2:1:2;
for k = 1:length(U)
    array = U{k};
    spikes = squeeze(array.R_ntk(1,:,:));
    touchIdx = [find(array.S_ctk(9,:,:)==1) ;find(array.S_ctk(12,:,:)==1)];
%     touchIdx = [find(array.S_ctk(9,:,:)==1)];
    blWin = touchIdx + swin; 
    postTpeak = touchIdx + peaks(k); 
    peakWin = postTpeak+swin; 
    
    bl = sum(spikes(blWin),2); 
    t = sum(spikes(peakWin),2); 

    
    change = t-bl; 
    
    figure(329);subplot(4,8,k)
    histogram(change,[-5:1:5])
    propTresp(k) = (sum(change>0)./numel(change)) - (sum(change<0)./numel(change));
end
    
    
    figure(490);clf
    histogram(propTresp,0:.05:1,'normalization','probability')
    set(gca,'ylim',[0 .5],'ytick',0:.25:.5,'xtick',0:.25:1)
    xlabel('proportion of touches responsive') 
    ylabel('number of cells') 