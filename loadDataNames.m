function [layer, txt] = loadDataNames(array)
U{1} = array{1};
if strcmp(U{1}.meta.layer,'L5b')
    [ndata, txt, alldata] =xlsread('CellsConversionChart170224','I23:J67');
    disp(U{1}.meta.layer);
    layer = 'L5b';
elseif strcmp(U{1}.meta.layer,'NL5b')
    [ndata, txt, alldata] =xlsread('CellsConversionChart170224','AB75:AC101');
    disp(U{1}.meta.layer);
    layer = 'NL5b';
elseif strcmp(U{1}.meta.layer,'L3')
    [ndata, txt, alldata] =xlsread('CellsConversionChart170224','B25:C44');
    disp(U{1}.meta.layer);
    layer = 'L3';
elseif strcmp(U{1}.meta.layer,'L4')
    [ndata, txt, alldata] =xlsread('CellsConversionChart170224','O23:P28');
    disp(U{1}.meta.layer);
    layer = 'L4';
elseif strcmp(U{1}.meta.layer,'L3Out')
    [ndata, txt, alldata] =xlsread('CellsConversionChart170224','B77:C82');
    disp(U{1}.meta.layer);
    layer = 'L3Out';
elseif strcmp(U{1}.meta.layer,'L5bOut')
    [ndata, txt, alldata] =xlsread('CellsConversionChart170224','I75:J82');
    disp(U{1}.meta.layer);
    layer = 'L5bOut';
else strcmp(U{1}.meta.layer,'L5bInt')
    [ndata, txt, alldata] =xlsread('CellsConversionChart170224','V75:W81');
    disp(U{1}.meta.layer);
    layer = 'L5bInt';
end
txt=txt(~isnan(ndata));
end