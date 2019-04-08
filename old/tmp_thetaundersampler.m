samp=3


jinho=wl.deltaKappa{1};
shortjinho=(jinho(1:samp:length(jinho)));
newjinho=zeros(length(shortjinho),2);
newjinho(:,1)=shortjinho;
newjinho(:,2)=1:samp:length(jinho);
figure;plot(jinho)
hold on
plot(newjinho(:,2),newjinho(:,1),'r')