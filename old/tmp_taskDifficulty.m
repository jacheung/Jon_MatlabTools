clear V
clear F1G RMSEgroup F1Gtree RMSEdecomp
[V] = classifierWrapper(U);

% PARAMETERS SETTING

clearvars -except V U BV D SM F1G RMSEgroup RMSEdecomp F1Gtree R PAS POP
numIterations = 10;

designvars = 'ubered';
% 1) 'theta' 2) 'pas' (phase amp midpoint) 3) 'counts' $) 'ubered'
classes = 'gonogo';
% 1) 'gonogo' 2) 'FAvsCR' 3) 'lick' 4) allBehavTypes

sample ='bias';
% 1) 'bias' (takes 70% from each class for train) 2) 'random' just takes
% random 70% to train

% Only for 'ubered' or 'pas'
normalization = 'none';
% 1) 'whiten' 2) 'none';

% Only for 'ubered'
removal = 'no';
% 1) 'yes' to remove 0 touch trials and auto classify as CR

balance = 'off';

taskD = zeros(1,length(U));
for rec = 1:length(U)
[DmatX, DmatY , motorX] = designMatrixBuilder(V(rec),U{rec},designvars,classes,normalization,removal,balance,'random');
      
gos = intersect(find(DmatY==1),find(DmatX(:,2)>0));
nogos = intersect(find(DmatY==2),find(DmatX(:,2)>0));

gothetasraw = DmatX(gos,1);
meanplusSD = mean(gothetasraw)+2*std(gothetasraw); %meanPLUS2SDs
gothetas=gothetasraw(gothetasraw<meanplusSD);

nogothetasraw = DmatX(nogos,1);
meanplusSDnogo = mean(nogothetasraw)-2*std(nogothetasraw);
nogothetas = nogothetasraw(nogothetasraw>=meanplusSDnogo);

% figure(2304);clf
% bar([-20:1:50],histc(gothetasraw,[-20:1:50]),'b')
% hold on;
% bar([-20:1:50],histc(nogothetasraw,[-20:1:50]),'r')

taskD(rec) = min(nogothetasraw) - max(gothetasraw);

[min(nogothetas) max(gothetas)]

    if strcmp(U{rec}.meta.layer,'D')
       POP.taskD{1}(rec,1) = min(nogothetas) - max(gothetas);
    elseif strcmp(U{rec}.meta.layer,'SM')
       POP.taskD{2}(rec,1) = min(nogothetas) - max(gothetas);
    elseif strcmp(U{rec}.meta.layer,'BV')
        POP.taskD{3}(rec,1) = min(nogothetas) - max(gothetas);
    end

end
