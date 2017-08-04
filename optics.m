[RD,CD,order]=tmp_optics([atTouchVar(:,1),atTouchVar(:,2)],10);
figure(5);bar(order,RD)

grab=order(414:440);

figure(6);clf;hold on
for i = 1:length(grab)
    scatter(atTouchVar(grab(i),1),atTouchVar(grab(i),2));
    hold on;
end
figure(7);clf;
spksper=sum(spikesAligned(grab,25:50),2);%calculating spks/touch for each point in chosen cluster

