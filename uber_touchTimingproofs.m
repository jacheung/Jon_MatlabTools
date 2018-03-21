[P] = findMaxMinProtraction(array,'avail2lick');
thetas = squeeze(array.S_ctk(1,:,:));

diffs = zeros(length(P.peakidx),1);
for i = 1:length(diffs)
    tmp = find(P.troughIdx<P.peakidx(i),1,'last');
    if ~isempty(tmp)
        diffs(i) = P.troughIdx(tmp);
    else
        diffs(i) = nan;
    end
end

match = [diffs P.peakidx];
[~, keep] = unique(match(:,1));
match = match(keep,:);


ttopeak= match(:,2)-match(:,1);
ttopeakthetas = thetas(match(:,2)) - thetas(match(:,1));

tosstrials = find(ttopeak>100);

ttopeak(tosstrials)=[];
ttopeakthetas(tosstrials)=[];

cutoffs = median(ttopeak)+2*std(ttopeak);
tossagain = find(ttopeak>cutoffs);
ttopeak(tossagain)=[];
ttopeakthetas(tossagain)=[];

troughwhisks = thetas(match(:,1));
figure(5420);clf;histogram(troughwhisks,'binedges',-20:1:20,'normalization','probability')
xlabel('Trough of whisk (theta)')
figure(234);clf;histogram(ttopeak,'binedges',0:1:cutoffs,'normalization','probability')
xlabel('Time from trough to peak whisk (ms)')

tperms= ttopeakthetas./ttopeak;
figure(20);clf;histogram(tperms,'binedges',0:.1:max(tperms),'normalization','probability')
xlabel('Theta traveled per ms during whisk')
%% Ttest plots for identifying touch and protraction start
thetas = squeeze(array.S_ctk(1,:,:));

xmin=1;
xmax=array.k;
n=1;
trial=round(xmin+rand(1,n)*(xmax-xmin));
trial =150
ranges = (trial*4000)+1:(trial*4000)+4000;

figure(58);clf;
plot(thetas(ranges(1):ranges(end)));

plotIdx = find(match(:,2)<ranges(end) & match(:,2)>ranges(1));

plotpeak=match(plotIdx,2);
plotonset=match(plotIdx,1);
hold on; scatter(plotpeak-ranges(1),thetas(plotpeak),'r')
hold on; scatter(plotonset-ranges(1),thetas(plotonset),'g')





