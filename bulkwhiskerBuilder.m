cd ('Z:\Users\Jon\DATA\WhiskerArrayBuilders\ToBuild')
filelist = dir('*.mat');

for i = 1:length(filelist)
    cd ('Z:\Users\Jon\DATA\WhiskerArrayBuilders\ToBuild')
    clearvars -except filelist i
    load(filelist(i).name);
    cd(d);
    
    Whisker.makeAllDirectory_WhiskerTrial(d,0,'mask', maskPoints,...
        'trial_nums',trialNums,'include_files',includef,...
        'barRadius',10,'faceSideInImage', 'top', 'framePeriodInSec',.001,...
        'imagePixelDimsXY',videoDim,'pxPerMm',33,'mouseName',mouseName,'sessionName',sessionName,'protractionDirection','leftward')
    
    Whisker.makeAllDirectory_WhiskerSignalTrial(d,'include_files',includef,'polyRoiInPix',[132-33 132+33],'follicleExtrapDistInPix',33);
    Whisker.makeAllDirectory_WhiskerTrialLiteI(d,'include_files',includef,'r_in_mm',4,'calc_forces',true,'whisker_radius_at_base', 36.5,'whisker_length', 18,'baseline_time_or_kappa_value',0);
    wl = Whisker.WhiskerTrialLiteArray(d);
    save([d mouseName sessionName '-WTLIA.mat'],'wl');
    
    tid = 0; % Set trajectory ID to view
    Whisker.view_WhiskerTrialLiteArray(wl,tid)
    delete(filelist(i).name)
end

for d = 1:length(filelist)
disp(['completed whiskerarrays for ' filelist(d).name ])
end