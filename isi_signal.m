function isi_signal(filename, chunksize,nchunks, baseline, stimdur)
% analyze the signal from isi image

nfiles = nargin;

if nfiles < 1
    error('Must input at least one file name.')
end

fn = filename;
signaloutput_fn = strcat(fn,'_signal.mat');
basePeriod = 1 : baseline;
stimPeriod = baseline + 1 : baseline + stimdur; 

load(strcat(fn,'_result.mat')) % need to be fixed with incremental result files...
diffMean = diffMean';
ylim_param = [-0.3 0.1];

figure, imagesc(diffMean, [-0.5 0.5]), axis image
radius = 30;
[x, y] = ginput(1);

mask = zeros(size(diffMean,1), size(diffMean,2));
for k = 1 : size(diffMean,1)
    for m = 1 : size(diffMean,2)
        if (m-x)^2 + (k-y)^2 < radius^2
            mask(k,m) = 1;
        end
    end
end

figure, imagesc(mask), axis image


f = 1;
tsimage_norm = zeros(size(diffMean,2), size(diffMean,1), chunksize, nchunks); %
for n = 1:nchunks
    rep = read_qcamraw([fn '.qcamraw'], f:(f+chunksize-1));
    rep = double(rep);
    temp_base = mean(rep(:,:,basePeriod),3);
    for p = 1 : chunksize
       tsimage_norm(:,:,p,n) = (rep(:,:,p) - temp_base)./temp_base * 100;
    end
    f = f+chunksize;
end


figure, subplot(1,3,1), imagesc(mean(mean(tsimage_norm(:,:,stimPeriod,1:5),4),3)', [-0.5 0.5]), axis image, title('1:5')
subplot(1,3,2), imagesc(mean(mean(tsimage_norm(:,:,stimPeriod,1:10),4),3)', [-0.5 0.5]), axis image, title('1:10')
subplot(1,3,3), imagesc(mean(mean(tsimage_norm(:,:,stimPeriod,1:20),4),3)', [-0.5 0.5]), axis image, title('1:20')

im_response = mean(mean(tsimage_norm(:,:,stimPeriod,:),4),3);

ts_all = zeros(nchunks,chunksize);
for q = 1 : nchunks
    for r = 1 : chunksize
        ts_all(q,r) = sum(sum(tsimage_norm(:,:,r,q)'.*mask))/sum(sum(mask));
    end
end
figure, subplot(1,3,1), plot(mean(ts_all(1:5,:),1)), ylim(ylim_param), title('1:5')
subplot(1,3,2), plot(mean(ts_all(1:10,:),1)),ylim(ylim_param), title('1:10')
subplot(1,3,3), plot(mean(ts_all(1:20,:),1)), ylim(ylim_param), title('1:20')

save(signaloutput_fn, 'tsimage_norm', 'ts_all', 'im_response')
