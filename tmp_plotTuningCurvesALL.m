    figure(44);clf;subplot(2,2,1)
    cdfplot(semi.thetaALL(:,1));
    hold on;f=cdfplot(cont.thetaALL(:,1));
    set(f,'Color','Cyan')
    set(gca,'xlim',[-30 60],'ytick',[0:.25:1])
    xlabel('Theta');ylabel('Percentile');title('');grid off
    
   
    subplot(2,2,2)
    cdfplot(semi.ampALL(:,1))
    hold on;f=cdfplot(cont.ampALL);
    set(f,'Color','Cyan')
    set(gca,'xlim',[0 50],'ytick',[0:.25:1])
        xlabel('Amplitude');ylabel('Percentile');title('');grid off
        hold on;legend('Semi n=11','Cont n=6')
        
    subplot(2,2,3)
    cdfplot(semi.spALL(:,1));  
    hold on; f=cdfplot(cont.spALL);
    set(f,'Color','Cyan')
    set(gca,'xlim',[-30 60],'ytick',[0:.25:1])
    xlabel('Setpoint');ylabel('Percentile');title('');grid off
    
    subplot(2,2,4)
    hold on;cdfplot(semi.phaseALL(:,1));
    hold on; f=cdfplot(cont.phaseALL);
    set(f,'Color','Cyan')
    set(gca,'xlim',[-pi pi],'xtick',pi*[-1:.5:1],'xticklabel',{'-pi','-pi/2',0,'pi/2','pi'},'ytick',[0:.25:1])
     xlabel('Phase');ylabel('Percentile');title('');grid off
     
%%
%REST vs NORM
range = [-30:50];
test=(histc(semi.thetaALL,range))./numel(semi.thetaALL);
testC=(histc(cont.thetaALL,range))./numel(cont.thetaALL);

figure(22);clf;subplot(2,1,1)
bar(range,test)
hold on;bar(range,testC,'c');alpha(.5)
hold on;legend('Semi n=11','Cont n=6')
xlabel('Theta')
ylabel('Proportion of Time')

test=(histc(semi.thetaNORM,range))./numel(semi.thetaNORM);
testC=(histc(cont.thetaNORM,range))./numel(cont.thetaNORM);

subplot(2,1,2);
bar(range,test)
hold on;bar(range,testC,'c');alpha(.5)
hold on;legend('Semi n=11','Cont n=6')
xlabel('Theta with Rest Shift')
ylabel('Proportion of Time')
%%
range = [-30:50];
test=(histc(semi.thetaALL,range))./numel(semi.thetaALL);
testC=(histc(cont.thetaALL,range))./numel(cont.thetaALL);

figure(22);clf;
subplot(2,2,1);bar(range,test)
hold on;bar(range,testC,'c');alpha(.5)
hold on;legend('Semi n=11','Cont n=6')
xlabel('Theta')
ylabel('Proportion of Time')

range = [0:50];
test=(histc(semi.ampALL,range))./numel(semi.ampALL);
testC=(histc(cont.ampALL,range))./numel(cont.ampALL);

subplot(2,2,2);bar(range,test)
hold on;bar(range,testC,'c');alpha(.5)
hold on;legend('Semi n=11','Cont n=6')
xlabel('Amplitude')
ylabel('Proportion of Time')

range = [-30:50];
test=(histc(semi.spALL,range))./numel(semi.spALL);
testC=(histc(cont.spALL,range))./numel(cont.spALL);

subplot(2,2,3);bar(range,test)
hold on;bar(range,testC,'c');alpha(.5)
hold on;legend('Semi n=11','Cont n=6')
xlabel('Setpoint')
ylabel('Proportion of Time')

range = [-3.14:.25:3.14];
test=(histc(semi.phaseALL,range))./numel(semi.phaseALL);
testC=(histc(cont.phaseALL,range))./numel(cont.phaseALL);

subplot(2,2,4);bar(range,test)
hold on;bar(range,testC,'c');alpha(.5)
hold on;legend('Semi n=11','Cont n=6')
xlabel('Phase')
ylabel('Proportion of Time')
