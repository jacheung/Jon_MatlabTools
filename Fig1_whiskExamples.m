mouseName = 'AH0712'
sessionName ='170804'
videoloc = 'JON'

d= (['Z:\Data\Video\' videoloc filesep mouseName filesep sessionName filesep])
load(['Z:\Users\Jon\DATA\BehaviorArrays\solo_' mouseName '_' sessionName '.mat'])

cd(d)

trialnums ={'136','134'};

for i = 1:length(trialnums)
    trialnum = trialnums{i};
WSTName = ['Z:\Data\Video\' videoloc filesep mouseName filesep sessionName filesep mouseName 'x' sessionName '-' trialnum '_' 'WST' '.mat'];
load(WSTName)%load file based on trial above to test mask 
tp = [.5 1.25];
%%%%%% plot any mask you want use trial number above  
figure(4329);subplot(1,2,i);ws.plot_fitted_whisker_time_projection(0,'k',tp)
hold on; 
ws.plot_fitted_whisker_ROI_time_projection(0,'r',tp)
ws.plot_mask(0,'g',tp);
ws.plot_follicle_position_time_projection(0,'b',tp)

tmp = load(['AH0712x170804-' trialnum '.bar']);
s=scatter(tmp(1,2),tmp(1,3),'filled','c');
s.SizeData = 200;

vid = VideoReader('AH0712x170804-4.mp4');

b.trials{str2double(trialnum)-1}.trialType
set(gca,'visible','off')
set(gca,'xlim',[0 vid.Width],'ylim',[0 vid.Height],'xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[])
end