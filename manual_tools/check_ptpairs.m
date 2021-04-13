clc;
addpath('align_config')
parentpath =  get_parent_path;
flipimg = true;
imgpath = [parentpath, filesep, 'stitched'];
secondpath = [imgpath, filesep, '2nd_round'];
imgpaths = {imgpath, secondpath};

kfpath = [parentpath, filesep, 'keyFrames'];
ptpairpath = [parentpath, filesep, 'ptpairs'];
maskrpath = [parentpath, filesep, 'masks_region'];

scl = 0.5;

imgptdir = [parentpath, filesep, 'imgpt'];
matlist = dir([ptpairpath, filesep, '*.mat']);
% matlist = dir([kfpath,filesep,'**',filesep,'ptpairs', filesep,'*.mat']);

% MATnames = erase({matlist.name}, '.mat');
% secondlist = dir([secondpath, filesep,'*.png']);
% idxt = contains(MATnames, erase({secondlist.name}, '.png'));
% idxt = count(MATnames, 'ref') < 2;
% idxt = contains({matlist.name}, 'W03_Sec028');
% matlist = matlist(idxt);

imglist = dir([imgpath, filesep, '*.png']);
imgnames_a = strrep({imglist.name},'.png','');

meshpath = [parentpath, filesep, 'meshes_affine'];
dis = nan(numel(matlist),1);
for k = 1:numel(matlist)
%     load([strrep(matlist(k).folder,'ptpairs0','ptpairs'), filesep, matlist(k).name]);
    load([matlist(k).folder, filesep, matlist(k).name]);
    yx1 = ptpairs.yx1;
    yx2 = ptpairs.yx2;
% 
%     M1 = load([meshpath, filesep, strrep(imgnames{1},'.png','.mat')], 'M');
%     M2 = load([meshpath, filesep, strrep(imgnames{2},'.png','.mat')], 'M');
%     M1 = M1.M;
%     M2 = M2.M;
%     N1 = size(M1.TR0.Points,1);
%     N2 = size(M2.TR0.Points,1);
%     A1 = geometries.fit_affine(fliplr(M1.TR.Points(1:10:N1,:)), fliplr(M1.TR0.Points(1:10:N1,:)));
%     A2 = geometries.fit_affine(fliplr(M2.TR.Points(1:10:N2,:)), fliplr(M2.TR0.Points(1:10:N2,:)));
%     yx1 = yx1 * A1(1:2,1:2) + A1(3,1:2);
%     yx2 = yx2 * A2(1:2,1:2) + A2(3,1:2);
% 
    A = geometries.fit_affine(yx1, yx2);
%     [U,S,V] = svd(A(1:2,1:2));
%     R = U * V';
    yx2 = yx2 * A(1:2,1:2) + A(3,1:2);
    dis0 = sum((yx2 - yx1).^2, 2).^0.5;
%     dis(k) = max(abs(log(diag(S))));
%     dis(k) = abs(R(2));
%     dis(k) = max(dis0);%/(std(dis0(:))+0.1);
%     dis(k) = range(maxk(dis0,5));
      dis(k) = -numel(yx1);
%     dis(k) = matlist(k).datenum;
%     dis(k) = -mean(ptpairs.conf);
end
[~, idx] = sort(dis, 'descend');
matlist = matlist(idx);

matnames = {matlist.name};
% 
% sttidx = find(strcmpi(matnames, 'sec2552_sec2553ref.mat'));
% endidx = find(strcmpi(matnames, 'sec2239ref_sec2240ref.mat'));
sttidx = 1;
endidx = numel(matlist);
load([matlist(sttidx).folder, filesep, matlist(sttidx).name]);
% IMG1 = imread([imgpath, filesep, imgnames{1}]);
% IMG2 = imread([imgpath, filesep, imgnames{2}]);
IMG1 = utils.try_read_img(imgpaths, imgnames{1});
IMG2 = utils.try_read_img(imgpaths, imgnames{2});
if scl ~= 1
    IMG1 = imresize(IMG1, scl);
    IMG2 = imresize(IMG2, scl);
end

cropbox2 = utils. draw_cropbox(IMG2);
if ~isempty(cropbox2)
    cropbox1 = utils.draw_cropbox(IMG1);
    if isempty(cropbox1)
        cropbox1 = cropbox2;
    end
else
    cropbox1 = [];
end

if exist([maskrpath, filesep, imgnames{1}], 'file')
    maskr1 = imread([maskrpath, filesep, imgnames{1}]);
    if scl ~= 1
        maskr1 = imresize(maskr1, scl, 'nearest');
    end
else
    maskr1 = ones(size(IMG1));
end
if exist([maskrpath, filesep, imgnames{2}], 'file')
    maskr2 = imread([maskrpath, filesep, imgnames{2}]);
    if scl ~= 1
        maskr2 = imresize(maskr2, scl, 'nearest');
    end
else
    maskr2 = ones(size(IMG2));
end

IMGs = {IMG1, IMG2};
maskrs = {maskr1, maskr2};
cropboxes = {cropbox1, cropbox2};

dyx1 = [0, 0];
dyx2 = [0, 0];

imgnames0 = imgnames;

tic;
for k = sttidx:1:endidx
    ppname = strrep(matlist(k).name,'.mat','');
    ppname = split(ppname,'_');
    tmpN = numel(ppname);
    sec1name = strjoin(ppname(1:(tmpN/2)), '_');
    sec2name = strjoin(ppname((tmpN/2+1):end), '_');
    ppname = {sec1name, sec2name};
    [~,ppidx] = ismember(ppname, imgnames_a);
    if abs(diff(ppidx)) > 1
        % continue
    end
    % if ~contains(matlist(k).name, 'Sec042')
        % continue;
    % end
    load([matlist(k).folder, filesep, matlist(k).name]);
    [~, idxt] = ismember(imgnames, imgnames0);
    if idxt(1) == 0
        IMG1 = utils.try_read_img(imgpaths, imgnames{1});
        if scl ~= 1
            IMG1 = imresize(IMG1, scl);
        end
        if exist([maskrpath, filesep, imgnames{1}], 'file')
            maskr1 = imread([maskrpath, filesep, imgnames{1}]);
            if scl ~= 1
                maskr1 = imresize(maskr1, scl, 'nearest');
            end
        else
            maskr1 = ones(size(IMG1));
        end
        if ~isempty(cropbox1)
            cropbox1 = utils.draw_cropbox(IMG1);
            if isempty(cropbox1)
                cropbox1 = cropbox2;
            end
        end
    else
        IMG1 = IMGs{idxt(1)};
        maskr1 = maskrs{idxt(1)};
        cropbox1 = cropboxes{idxt(1)};
    end
    if idxt(2) == 0
        IMG2 = utils.try_read_img(imgpaths, imgnames{2});
        if scl ~= 1
            IMG2 = imresize(IMG2, scl);
        end
        if exist([maskrpath, filesep, imgnames{2}], 'file')
            maskr2 = imread([maskrpath, filesep, imgnames{2}]);
            if scl ~= 1
                maskr2 = imresize(maskr2, scl, 'nearest');
            end
        else
            maskr2 = ones(size(IMG2));
        end
        if ~isempty(cropbox2)
            cropbox2 = utils.draw_cropbox(IMG2);
            if isempty(cropbox2)
                cropbox2 = cropbox1;
            end
        end
    else
        IMG2 = IMGs{idxt(2)};
        maskr2 = maskrs{idxt(2)};
        cropbox2 = cropboxes{idxt(2)};
    end
    IMGs = {IMG1, IMG2};
    maskrs = {maskr1, maskr2};
    cropboxes = {cropbox1, cropbox2};
    imgnames0 = imgnames;
    
    if ~isempty(cropbox1)
        y0 = max(1, cropbox1(1));
        y1 = cropbox1(2);
        x0 = max(1, cropbox1(3));
        x1 = max(1, cropbox1(4));
        IMG1t = IMG1(y0:min(y1,size(IMG1,1)),x0:min(x1,size(IMG1,2)));
        maskr1t = maskr1(y0:min(y1,size(maskr1,1)),x0:min(x1,size(maskr1,2)));
        dyx1 = [y0, x0] - 1;

        y0 = max(1, cropbox2(1));
        y1 = cropbox2(2);
        x0 = max(1, cropbox2(3));
        x1 = max(1, cropbox2(4));
        
        IMG2t = IMG2(y0:min(y1,size(IMG2,1)),x0:min(x1,size(IMG2,2)));
        maskr2t = maskr2(y0:min(y1,size(maskr2,1)),x0:min(x1,size(maskr2,2)));
        dyx2 = [y0, x0] - 1;
        if all(maskr1t(:) == 0) || all(maskr2t(:) == 0)
            continue
        end
    else
        IMG1t = IMG1;
        IMG2t = IMG2;
        maskr1t = maskr1;
        maskr2t = maskr2;
        dyx1 = [0, 0];
        dyx2 = [0, 0];
    end
    disp([num2str(k), ': ',matlist(k).name]);
    ptpairt = ptpairs;
    ptpairt.yx1 = ptpairt.yx1 * scl - dyx1;
    ptpairt.yx2 = ptpairt.yx2 * scl - dyx2;

    if flipimg
        [ptpairt, modified] = edit_ptpairs(IMG2t, IMG1t, geometries.flip_ptpairs(ptpairt), {maskr2t, maskr1t});
        ptpairt = geometries.flip_ptpairs(ptpairt);
    else
        [ptpairt, modified] = edit_ptpairs(IMG1t, IMG2t, ptpairt, {maskr1t, maskr2t});
    end
    if modified == -1
        disp('exit manual ptpair check.')
        return
    end
    if modified == -2
        delete([matlist(k).folder, filesep, matlist(k).name])
        if contains(ptpairpath, 'ptpairs0')
            try
                delete([strrep(matlist(k).folder, 'ptpairs0', 'ptpairs'), filesep, matlist(k).name]);
                disp(['deleted ', matlist(k).name])
            end
        end
        continue
    end
    ptpairt.yx1 = ptpairt.yx1 + dyx1;
    ptpairt.yx2 = ptpairt.yx2 + dyx2;
    ptpairs = ptpairt;
    ptpairs.yx1 = ptpairs.yx1/scl;
    ptpairs.yx2 = ptpairs.yx2/scl;
    clf
    if modified
        save([matlist(k).folder, filesep, matlist(k).name], 'ptpairs','imgnames');
        if contains(ptpairpath, 'ptpairs0')
            try
                delete([strrep(matlist(k).folder, 'ptpairs0', 'ptpairs'), filesep, matlist(k).name]);
                disp(['deleted ', matlist(k).name])
            end
        end
    end
    % return
    % saveas(gcf,[imgptdir, filesep, strrep(matlist(k).name,'.mat','.png')])
end

close all