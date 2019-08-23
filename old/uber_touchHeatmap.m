%%Touch Heatmap
totonset=cell(1,length(L3)); %each cell = 1 neuron with spks/touch. Timepoints from -50:150ms after touch
for rec=1:length(L3)
    touchIdx = [find(L3{rec}.S_ctk(9,:,:)==1);find(L3{rec}.S_ctk(12,:,:)==1)];
    spikesAligned = zeros(numel(touchIdx),75);
    spikes = squeeze(L3{rec}.R_ntk);
 
    for i = 1:size(spikesAligned,1)
        spikesAligned(i,:) = spikes(touchIdx(i)+[-25:49]);
    end

totonset{rec}=(sum(spikesAligned)/numel(touchIdx))';
end

%%normalize
totnorm=cell(1,length(totonset));
for i=1:length(totonset)
    small=min(totonset{i});
maxmin=max(totonset{i})-min(totonset{i});
totnorm{i}=(totonset{i}-small)/maxmin;
end
figure;
imagesc(cell2mat(totnorm)')
hold on
plot([25 25],[25 0],'w:') 
title('L3 Touch Onset')
set(gca,'YTickLabel',[])
set(gca,'YTick',[])
set(gca,'XTick',[])
%set(gca,'XTickLabel',{'-25' ;'0' ; '25' ; '50'})
hCBar=colorbar;
set(hCBar,'yaxislocation','right')
set(hCBar,'YTick',[1:1:10])