function isi_register(fn, varargin)

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
thresh = 0.1; % default threshold
hsize = 10;
sigma = 3;
h = fspecial('gaussian', hsize, sigma);
radius = 20; % defaul radius
%  S = regionprops(CC, 'Area');
% [maxarea, index] = max(S); % automatically selecting the component with the larges area. subject to modification...

%%


vas_fn = 'vas.qcamraw';
def_fn = 'def.qcamraw';
if ~exist(vas_fn)
    error('''vas.qcamraw'' file is required. Image the vasculature');
elseif ~exist(def_fn)
    warning('''def.qcamraw'' is not found. Proceeding with red image...');
    def_fn = strcat(fn,'.qcamraw');
end
vas = read_qcamraw(vas_fn,1);
def = read_qcamraw(def_fn,1);
% vas = imread('vas.tif');
% def = imread('def.tif');


mode_set = ['auto', 'manual', 'circ'];
switch (nargin)
    case 1,
        mode = 'auto';
    case 2,
        if isnumeric(varargin{1})
            if varargin{1} < 1 && varargin{1} > 0
                mode = 'auto';
                thresh = varargin{1};
            else
                error('threshold should be between 0 and 1')
                return
            end
        else
            mode = varargin{1};
            if isempty(strfind(mode_set,mode))
                error('mode should be ''auto'', ''manual'', or ''circ''');
                return
            end
        end
    case 3,
        mode = varargin{1};
        if isempty(strfind(mode_set,mode))
            error('mode should be ''auto'', ''manual'', or ''circ''');
        end
        
        if strcmp(mode,'circ')
            if ~isnumeric(varargin{2})
                error('3rd argument should be an number');
            else
                radius = round(varargin{2});
            end
        else
            if varargin{2} < 1 && varargin{2} > 0
                thresh = varargin{2};
            else
                error('threshold should be between 0 and 1')
                return
            end
        end
    otherwise,
        error('too many input arguments!!');
end

load_fn = strcat(fn, '_result.mat'); 
if ~exist(load_fn)
    error('''_result.mat'' file is required. Run ''isi_image()''');
else

load(load_fn, 'diffMean'); 

im_filt = imfilter(diffMean,h);

%
%
%
% hsize and sigma for guassian filter is set to 3 and 0.5. 
% hsize = 3;
% sigma = 0.5;
%
%
if strcmp(mode,'auto')
    %
    %
    % manual selection of the threshold
%     thresh = 0.1;
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

elseif strcmp(mode,'manual')
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
else % strcamp(mode,'circ')
    figure, imagesc(im_filt'), axis image;
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

% calculating registration function
[optimizer, metric] = imregconfig('monomodal');
tform = imregtform(def,vas,'rigid',optimizer,metric);

% registration using the transform function
im_roi = imwarp(im_roi,tform);
 
tp_vas = vas';
tp_roi = bwperim(im_roi');
tp_roi = imresize(tp_roi,[size(vas,2) size(vas,1)]);
% showing the result
im_fused = imfuse(tp_vas,tp_roi,'diff');
imtool(im_fused);

save_fn = strcat(fn,'_register.mat');
save(save_fn, 'im_fused','mode','thresh','radius', 'hsize', 'sigma','vas','def');
if exist('im_filt')
    save(save_fn, 'im_filt', '-append');
end
imwrite(im_fused, 'registered.tif');
end


