function isi_overlapOLD(fn, vsfn, varargin)

% Calculates intensity-based image transformation (rigid) from fn_def 
% (filename of defocused vessel image, usually 'def.qcamraw') to fn_vas 
% (filename of vasculature image, usually 'vas.qcamraw'), and then apply
% that transformation function to register roi to vasculature image
%
% Requires 'fn_signal.mat' file calculated by "isi_signal.m" which contains
% roi image (tsimage_norm)
% 
% mode can be 'auto', 'manual', and 'circ'. 'auto'-mode draws roi and 
% picks up automatically: good for when there is only one roi popped up, 
% 'manual'-mode requires selection among many rois popped up or when the
% roi does not seem to be reasonable, and 'circ' requires manual selection
% of the roi centroid and radius of a circle - default is 20 (pixels)
% Default mode is 'auto'
%
% JKim 01/10/2016

%% some variables manually set for now...
thresh = 0.01;
hsize = 10;
sigma = 3;
radius = 10;
%  S = regionprops(CC, 'Area');
% [maxarea, index] = max(S); % automatically selecting the component with the larges area. subject to modification...

%%

mode_set = ['auto', 'manual', 'circ'];
switch (nargin)
    case 2,
        mode = 'auto';
    case 3,
        mode = varargin{1};
        if isempty(strfind(mode_set,mode))
            error('mode should be ''auto'', ''manual'', or ''circ''');
        end
    case 4,
        mode = varargin{1};
        if isempty(strfind(mode_set,mode))
            error('mode should be ''auto'', ''manual'', or ''circ''');
        end
        
        if strcmp(mode, 'circ')
%             if ~isinteger(varargin{2})
%                 error('4th argument should be an integer for the radius');
%             else
                radius = varargin{2};
%             end
        else
            thresh = varargin{2};
        end
    otherwise,
        error('too much input argument');
end

load_fn = strcat(fn, '_result.mat'); 
if ~exist(load_fn)
    error('''_result.mat'' file is required. Run ''isi_image()''');
else

load(load_fn, 'diffMean'); 
 vas_fn = strcat(vsfn,'.qcamraw');
 vas = read_qcamraw(vas_fn,1);
%vas_fn = strcat(vsfn,'.tif');
%vas = imread(vas_fn);
vas = vas(:,:,1);
vas = vas';
%
%
%
% hsize and sigma for guassian filter is set to 3 and 0.5. 
% hsize = 3;
% sigma = 0.5;
%
%
if strcmp(mode,'auto')
    h = fspecial('gaussian', hsize, sigma);
    im_filt = imfilter(diffMean,h);
    %
    %
    % manual selection of the threshold
%     thresh fn= 0.1;
    %
    %
    %
    im_bin_res = im_filt < -thresh;
    im_bin_res = imfill(im_bin_res,'holes');
    CC = bwconncomp(im_bin_res);
    S = regionprops(CC, 'Area');
    sa = cat(1, S.Area);
    [maxarea, index] = max(sa); % automatically selecting the component with the larges area. subject to modification...
    im_roi = zeros(size(diffMean));
    im_roi(CC.PixelIdxList{1, index}) = 1;
    im_roi = im_roi';

elseif strcmp(mode,'manual')
    h = fspecial('gaussian', hsize, sigma);
    im_filt = imfilter(diffMean,h);
    im_bin_res = im_filt < -thresh;
    im_bin_res = imfill(im_bin_res,'holes');
    
    imshowpair(vas',im_filt','montage')
    figure, im_roi = bwselect(im_bin_res');
    im_roi = im_roi';
%     CC = bwconncomp(im_bin_res);
%     
%     figure, subplot(1,2,1), imagesc(im_filt), subplot(1,2,2), imagesc(im_bin_res)
%     while (1)
%        [x, y] = ginput(1); % select one point of the response image
%        if im_bin_res(x,y) == 0
%            warning('select the point inside one of the response regions');
%        else
%            exit
%        end
%     end
else % mode == 'circ'
    figure, imagesc(diffMean'), axis image;
    [x, y] = ginput(1);
    im_roi = zeros(size(diffMean,1), size(diffMean,2));
    for k = 1 : size(diffMean,1)
        for m = 1 : size(diffMean,2)
            if (k-x)^2 + (m-y)^2 < radius^2
                im_roi(k,m) = 1;
            end
        end
    end
end

% registration using the transform function
 
tp_vas = vas';
tp_roi = bwperim(im_roi');
tp_roi = imresize(tp_roi,[size(vas,2) size(vas,1)]);
% showing the result
im_fused = imfuse(tp_vas,tp_roi,'diff');
imtool(im_fused);

save_fn = strcat(fn,'_register.mat');
save(save_fn, 'im_fused','mode','thresh','radius', 'hsize', 'sigma','vas');
if exist('im_filt')
    save(save_fn, 'im_filt', '-append');
end
imwrite(im_fused, 'registered.tif');
 
end


