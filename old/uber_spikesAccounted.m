
for rec = 1:length(U)

spikes = squeeze(U{rec}.R_ntk);
touchPeriods = [find(U{rec}.S_ctk(11,:,:)==1) ;find(U{rec}.S_ctk(14,:,:)==1)];
touchPeriods = unique(touchPeriods);

whiskPeriods = [find(U{rec}.S_ctk(3,:,:)>=2.5)];
[~,~,removals] = intersect(touchPeriods,whiskPeriods);
whiskPeriods(removals) = [];

quietPeriods = [find(U{rec}.S_ctk(3,:,:)<2.5)];
[~,~,removals2] = intersect(touchPeriods,quietPeriods);
quietPeriods(removals2)=[];



totSpikeCount = sum(sum(spikes));

quietSpikeCount = sum(spikes(whiskPeriods));
whiskSpikeCount = sum(spikes(quietPeriods));
touchSpikeCount = sum(spikes(touchPeriods));


accounted = (sum([quietSpikeCount whiskSpikeCount touchSpikeCount]) / totSpikeCount)*100;

disp([ num2str(accounted) '% of spikes accounted for']);

figure(23);subplot(5,4,rec)
pie([quietSpikeCount whiskSpikeCount touchSpikeCount],{'Quiescent','Whisk','Touch'})

figure(25);subplot(5,4,rec)
pie([numel(quietPeriods) numel(whiskPeriods) numel(touchPeriods)],{'Quiescent','Whisk','Touch'})
end
% 
% figure(4);clf
% 
% for i = 1
%     
% touchPeriods = [find(U{rec}.S_ctk(11,:,i)==1) ;find(U{rec}.S_ctk(14,:,i)==1)];
% touchPeriods = unique(touchPeriods);
%     whiskPeriods = [find(U{rec}.S_ctk(3,:,i)>=2.5)];
%     [~,~,removals] = intersect(touchPeriods,whiskPeriods);
%     whiskPeriods(removals) = [];
%     quietPeriods = [find(U{rec}.S_ctk(3,:,i)<2.5)];
%     [~,~,removals2] = intersect(touchPeriods,quietPeriods);
%     quietPeriods(removals2)=[];
%     
%     if sum(U{rec}.R_ntk(1,:,i))>0
%         plot(find(U{rec}.R_ntk(1,:,i)==1),i,'k.')
%         hold on; plot(touchPeriods,i,'y.')
%         hold on; plot(whiskPeriods,i,'b.')
%         hold on; plot(quietPeriods,i,'c.')
%     else
%          plot(touchPeriods,i,'y.')
%         hold on; plot(whiskPeriods,i,'b.')
%         hold on; plot(quietPeriods,i,'c.')
% 
%     end
%     
% end
% set(gca,'ylim',[0 U{rec}.k]+1)