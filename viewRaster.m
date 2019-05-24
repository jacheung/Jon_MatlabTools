function viewRaster(array)

if nargin<1
    error('input one U array candidate')
end

%optional raster of FRs for tuned cells. 
figure(88);clf
allSpks = squeeze(array.R_ntk);
[~,idx] = sort(array.meta.motorPosition);
allSpks = allSpks(:,idx);
for k = 1:size(allSpks,2)
    st = find(allSpks(:,k)==1);
    if ~isempty(st)
        figure(88);hold on
        scatter(st,ones(length(st),1).*k,[],'.k')
    end
end
set(gca,'ydir','reverse','ylim',[1 array.k],'xtick',0:1000:4000)
ylabel('sorted trials from far to near')
xlabel('time from trial start (ms)')

pOnset = round(mean(array.meta.poleOnset)*1000);
hold on; plot([pOnset pOnset],[0 array.k],'-.g')

% spOffset = round(array.meta.poleOffset(idx)*1000); 
% spOffset(spOffset>array.t) = array.t;
% hold on; scatter(spOffset,1:size(allSpks,2),'ro')
% 




