function [norm,avg,range] = featureNormal(h)
    
avg = mean(h);
range = (max(h)-min(h));
range = repmat(range,size(h,1),1);
norm = (h - mean(h)) ./ range;
