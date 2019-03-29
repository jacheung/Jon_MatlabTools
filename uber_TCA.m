Uraw=[NL5b BVL5b];

clearvars -except BVL5b NL5b Uraw
touchCells = touchCell(Uraw,2,.5);
selectedCells = find(touchCells==1);
%%
 U = Uraw(selectedCells);
% U = Uraw;

window= [-25:50]; %spike sum from 25ms after touch idx 
touchStack = nan(length(U),length(window),5000);

for rec = 1:length(U)
touchIdx = [find(U{rec}.S_ctk(9,:,:)==1);find(U{rec}.S_ctk(12,:,:)==1)];
spikes = squeeze(U{rec}.R_ntk);
spikesAtTouch = spikes(touchIdx+window);
thetas = squeeze(U{rec}.S_ctk(1,:,:));
amp = squeeze(U{rec}.S_ctk(3,:,:));
mp = squeeze(U{rec}.S_ctk(4,:,:));
phase = squeeze(U{rec}.S_ctk(5,:,:));

ttouch{rec} = [thetas(touchIdx) amp(touchIdx) mp(touchIdx) phase(touchIdx)] ; 

for k = 1:length(touchIdx)
    touchStack(rec,1:length(window),k) = smooth(spikesAtTouch(k,:),10);
end

end
        %% 
close all
% convert data to a tensor object
cellNum = 15;
chosenStack = touchStack(cellNum,:,:);
nants = find(~isnan(chosenStack(1,1,:)));
chosenStack = chosenStack(1,:,nants);

tdata = tensor(chosenStack);

% rank of the data. (the number of latent factors.)
R = 1:10;
% fit the cp decomposition from random initial guesses
n_fits = 30;
 err = zeros(length(R),n_fits);
plotVals = cell(3,1);
 
for Rvals = 1:length(R) 
for n = 1:n_fits
    % fit model
    est_factors = cp_als(tensor(tdata),R(Rvals));
    
    % store error
    err(Rvals,n) = norm(full(est_factors) - tdata)/norm(tdata);
    
%     visualize fit for first several fits
%     if n < 4
%         % score aligns the cp decompositions
% %         [sc, est_factors] = score(est_factors, true_factors);
% %         
%         % plot the estimated factors
%         [~,plotVals{n} ] =viz_ktensor(est_factors, ... 
%             'Plottype', {'bar', 'line', 'scatter'}, ...
%             'Modetitles', {'neurons', 'time', 'trials'});
%         set(gcf, 'Name', ['estimated factors - fit #' num2str(n)]);
%     end
end
end


realTfactors = normalize_var(ttouch{cellNum},0,1);




modelTfactors = cell(3,1);
for iteration = 1:length(info) 
    currFactors = nan(R,length(chosenStack));
    modelTtime{iteration} = plotVals{iteration}.time;
    modelTfactors{iteration} = normalize_var(plotVals{iteration}.trials,0,1);
    valCorr{iteration} = corr([modelTfactors{iteration} realTfactors]);
end

modelnum =3; 
figure(588);clf
for compnum = 1:4
subplot(2,4,compnum);
hold on; plot(window,modelTtime{modelnum}(:,compnum),'k')
set(gca,'xtick',-25:25:50,'xlim',[-25 50])

compdef = valCorr{modelnum};
subplot(2,4,compnum+4);
bar(1:4,compdef(5:end,compnum))
end



figure(34);clf;hold on
plot(0:Rvals,[1; mean(err,2)],'k')
% plot(randn(n_fits,1), err, 'ob')
set(gca,'ylim',[0 1],'ytick',0:.5:1,'xtick',0:5:10)
ylabel('error');xlabel('# of components')


