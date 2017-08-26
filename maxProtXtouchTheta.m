% inputs = {hx_prevT,mx_prevT,FAx_prevT,CRx_prevT};

    hx_prevT = find(V(rec).trialNums.matrix(1,:)==1); % using current T num b/c licks shifted by padded 0
    mx_prevT = find(V(rec).trialNums.matrix(2,:)==1);
    FAx_prevT = find(V(rec).trialNums.matrix(3,:)==1); % using current T num b/c licks shifted by padded 0
    CRx_prevT = find(V(rec).trialNums.matrix(4,:)==1); % using current T num b/c licks shifted by padded 0
    


test=[]
tmp = {'CR'};
inputs = {CRx_prevT}
thetas = V(rec).var.(tmp{1}){1};
Ttmp=[V(rec).touchNum.(tmp{1})];

for g = 1:length(tmp) %for trial types...
    for k = 1:length(inputs{g}) %for number of T within trial types...
        d = inputs{g}(k); %for specific trial number...
        test(k) = median(P.theta(P.trialNums==d));
    end
end
% test(k) = mean(P.theta(P.trialNums==d));
meanMaxProt = test;
meanMaxProt(meanMaxProt == -50) = NaN

%MEAN TOUCHES

yesTouches = find(Ttmp>0);
noTouches = find(Ttmp==0);
tmp2 = [Ttmp zeros(1,length(Ttmp))];
for i = 1:length(Ttmp)
    tmp0 = [zeros(1,i) Ttmp zeros(1,length(Ttmp)-i)];
    tmp2 = [tmp2 ;tmp0];
end
tIdx = sum(tmp2(:,1:length(Ttmp)));
tIdx = [1 unique(tIdx)];
meanTtheta = NaN(1,length(Ttmp));
for i = 1:length(tIdx)-1
    meanTtheta(yesTouches(i)) = median(thetas(tIdx(i):tIdx(i+1)));
end

[meanMaxProt ; meanTtheta]

howfar=diff([meanMaxProt ; meanTtheta]);

figure(38);clf;scatter(meanTtheta,meanMaxProt,'b')
hold on; plot([-10 60],[-10 60],'-.k')
xlabel('Mean Theta at Touch/Trial')
ylabel('Max Protract')
axis('square')
hold on; text(25,0,['meandiff = ' num2str(nanmean(howfar))])