% S1datasetWrapper.m compiles an intermediate data structure (L) required
% for the matlab figure building scripts functions to properly operate.
%
% The script should be executed in a single directory where the individual
% recording data and meta_data files reside.
%
% SAHires 160714
cd('Z:\Data\ForDistribution\FINAL\matlab_scripts')
cd ..
cd datafiles

d = dir('data*.mat')
m = [16 14 27 45 48 17 32 28 39 20];
for i = 1:length(m)
    load(d(m(i)).name)
    L{i}.varnames = {    'thetaAtBase'    'unused'    'amplitude'    'setpoint'    'phase'    'deltaKappa'    'unused'    'unused'    'firstTouchOnset'    'firstTouchOffset'...
        'firstTouchAll'    'lateTouchOnset'    'lateTouchOffset'    'lateTouchAll'    'PoleAvailable'    'beamBreakTimes'};
    L{i}.cellNum = str2num(d(i).name(20:21));
    
    L{i}.t = sum(c.timeSeriesArrayHash.value{1}.trial == c.timeSeriesArrayHash.value{1}.trial(1));
    L{i}.k = length(c.trialIds);
    L{i}.c = 16;
    L{i}.u = 1;
    L{i}.S_ctk = zeros(L{i}.c, L{i}.t, L{i}.k);
    L{i}.S_ctk([9:14 16],:,:) = NaN;
    
    L{i}.R_ntk = zeros(1, L{i}.t, L{i}.k);
    
    Lmap = [1 3 5 4 0 6 0 0 0 15 16];
    for j = [1 2 3 4 6 10 11];
        L{i}.S_ctk(Lmap(j),:,:) = reshape(c.timeSeriesArrayHash.value{1}.valueMatrix(j,:),L{i}.t,L{i}.k) ;
    end
    
    Ton = reshape(c.timeSeriesArrayHash.value{1}.valueMatrix(7,:),L{i}.t,L{i}.k);
    Toff = reshape(c.timeSeriesArrayHash.value{1}.valueMatrix(8,:),L{i}.t,L{i}.k);
    
    for j = 1:size(Ton,2)
        L{i}.S_ctk(12,find(Ton(:,j)==1),j) = 1;
        L{i}.S_ctk(9,find(Ton(:,j)==1,1),j) = 1;
        L{i}.S_ctk(12,find(Ton(:,j)==1,1),j) = NaN;
        
    end
    
    for j = 1:size(Toff,2)
        L{i}.S_ctk(13,find(Toff(:,j)==1),j) = 1;
        L{i}.S_ctk(10,find(Toff(:,j)==1,1),j) = 1;
        L{i}.S_ctk(13,find(Toff(:,j)==1,1),j) = NaN;
        
    end
    
    tmp = nan(size(Ton));
    Ton_idx = find(Ton);
    Toff_idx = find(Toff);
    
    for j = 1:sum(Ton(:))
        
        tmp(Ton_idx(j):Toff_idx(j)) = 1;
    end
    
    L{i}.S_ctk(14,:, :) = tmp;
    
    
    for j = 1:size(Ton,2)
        
        L{i}.S_ctk(11,find(Ton(:,j)==1,1):find(Toff(:,j)==1,1),j) = 1;
        L{i}.S_ctk(14,find(Ton(:,j)==1,1):find(Toff(:,j)==1,1),j) = NaN;
    end
    
    % Account for 10ms whisker trial time offset in SAH recordings and different time format between SAH & JY recordings.
    if  strcmp(d(i).name(40:41),'JY')
        L{i}.R_ntk(1,:,:) = reshape(sum(reshape([c.timeSeriesArrayHash.value{2}.valueMatrix(2,6:end) 0 0 0 0 0],10,length(c.timeSeriesArrayHash.value{2}.valueMatrix(2,:))/10)), L{i}.t, L{i}.k);
    else
        L{i}.R_ntk(1,:,:) = reshape(sum(reshape([c.timeSeriesArrayHash.value{2}.valueMatrix(2,5:end)  0 0 0 0],10,length(c.timeSeriesArrayHash.value{2}.valueMatrix(2,:))/10)), L{i}.t, L{i}.k);
        
        L{i}.R_ntk(1,1:end-10,:) = L{i}.R_ntk(1,11:end,:);
        L{i}.R_ntk(1,end-9:end,:) =0;
               
    end
    
    %METADATA
    L{i}.meta.poleOnset = c.trialPropertiesHash.value{2}; %pole onset time
    L{i}.meta.poleOffset = c.trialPropertiesHash.value{3}; %pole offset time
    L{i}.meta.trialType = sum(c.trialTypeMat([1 2],:)); %trial type but summing all hit and miss = all go trials
    L{i}.meta.trialCorrect =sum(c.trialTypeMat([1 3],:)); %sum of all hit and CR trials = all trial correct
    L{i}.meta.layer = 'D';
    
    L{i}.trialTypeMat = c.trialTypeMat
    aux = findstr(d(m(i)).name,'_');
    strCell = d(m(i)).name(aux(2)+1:aux(3)-1);
    metaDataFileName = ['metadata_structure_' strCell '.mat'];
    load(metaDataFileName);
    L{i}.onsetLatency = meta_data.onsetLatency;
end
