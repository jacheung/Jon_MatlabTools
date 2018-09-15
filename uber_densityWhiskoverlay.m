%% load SEMI
load('Z:\Users\Jon\Projects\Characterization\SM.mat');
mouseNum = {'AH0349','AH0378','AH0407','AH0456','AH0454', 'AH0459','AH0374','AH0347','AH0281','AH0458'} ;
session = {'160511','160606','160613','160719','160706','160727','160309','160224','160212','160802'};

%% load BV
load('Z:\Users\Jon\Projects\Characterization\BV.mat');
mouseNum = {'AH0641','AH0656','AH0667','AH0638','AH0669', 'AH0657','AH0668','AH0717','AH0716','AH0712'} ;
session = {'170303','170228','170317','170214','170317','170314','170316','170901','170816','170804'};

%%
for rec = 2
    
    farthestnogoT= find(U{rec}.meta.motorPosition<min(U{rec}.meta.motorPosition)+10000);
    [objmask] = assist_touchmasks(U{rec});
    samplingIdx=find(objmask.samplingp(:,1)==1);
    ampthresh = squeeze(U{rec}.S_ctk(3,:,:)>5);
    cd(['Z:\Data\Video\JON\' mouseNum{rec} filesep session{rec} ])
    xcollat = []; ycollat = [];
    for k = 1:length(farthestnogoT)
        trialNum = num2str(farthestnogoT(k));
        if exist(['Z:\Data\Video\JON\' mouseNum{rec} filesep session{rec} filesep mouseNum{rec} 'x' session{rec} '-' trialNum '_WST.mat'],'file')
            load(['Z:\Data\Video\JON\' mouseNum{rec} filesep session{rec} filesep mouseNum{rec} 'x' session{rec} '-' trialNum '_WST.mat'])
            
            [x, y] = ws.plot_fitted_whisker_time_projection(0,'k');
            ampIdx = find(ampthresh(:,farthestnogoT(k))==1);
            selectedIdx = intersect(ampIdx,samplingIdx);
            
            if ~isempty(selectedIdx)
                xcollat = [xcollat [x{selectedIdx}]];
                ycollat = [ycollat [y{selectedIdx}]];
            end
        end
    end
end

densityplot(xcollat,ycollat,'nbins',[200 150])
cmap = colormap;
colormap(cmap(end:-1:1,:))
