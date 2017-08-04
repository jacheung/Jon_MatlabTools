function [range]=assist_bestBRange(T,span)
%insert trial array to output best range of span (# of trials you want
%calculated). Chooses the first best span of trials as there could be
%multiple ones. 

total=zeros(1,length(T.trialCorrects)-span);
for i = span:length(T.trialCorrects)
    total(i-(span-1))=sum(T.trialCorrects(i-(span-1):i));
end
best=find(total==max(total));
range=([best(1) best(1)+span]);
