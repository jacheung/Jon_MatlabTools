for rec=1:length(U) %for each cell
    [touchespredecision,GprelixtimesGO,GprelixtimesNOGO,Gprelixpcth,Gprexliduration,Gprelixtimes]...
        = assist_predecisionVar(U,rec);

wo=cell(1,length(Gprelixtimes)); %find whisker occupancy time BEFORE 1st touch
rest=cell(1,length(Gprelixtimes)); %find nonwhisking parts before first touch
for i=1:length(Gprelixtimes)
    if isempty(Gprelixtimes{i})
    whiskwin=find(U{rec}.S_ctk(3,:,i)>2.5);
    nowhiskwin=find(U{rec}.S_ctk(3,:,i)<2.5);
    wo{i}=U{rec}.S_ctk(1,whiskwin(1:end),i);
    rest{i}=U{rec}.S_ctk(1,nowhiskwin(1:end),i);
    else 
    whiskwin=find(U{rec}.S_ctk(3,1:Gprelixtimes{i}(1),i)>2.5);%determine amp>2.5 prefirst touch
    nowhiskwin=find(U{rec}.S_ctk(3,1:Gprelixtimes{i}(1),i)<2.5);
    wo{i}=U{rec}.S_ctk(1,whiskwin(1:end),i);
    rest{i}=U{rec}.S_ctk(1,nowhiskwin(1:end),i);
    end
end
gotrials=U{rec}.meta.trialType==1;
nogotrials=U{rec}.meta.trialType==0;

gowo=wo(gotrials);
gorest=rest(gotrials);
nogowo=wo(nogotrials);
nogorest=rest(nogotrials);

figure(34);clf;plot(-75:75,smooth(histc(cell2mat(gowo),-75:75),5))
hold on; plot(-75:75,smooth(histc(cell2mat(nogowo),-75:75),5),'r')
hold on; plot(-75:75,smooth(histc(cell2mat(rest),-75:75),5),'k')
xlabel('Whisker theta (deg)')
ylabel('Occupancy (ms)')

print(gcf,'-dpng',['Z:\Users\Jon\Projects\Characterization\' layer '\Figures\Strategy\' num2str(rec) '_' 'WhiskingOccupancy'])

end

%% Lick P vs Motor Pos 