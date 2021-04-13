clc;
addpath('align_config')
parentpath =  get_parent_path;
imgpath = [parentpath, filesep, 'stitched'];
discardpath = [imgpath, filesep, '2nd_pass'];
ptpairpath = [parentpath, filesep, 'ptpairs'];
outdir = [ptpairpath, filesep, 'covered'];
if ~exist(outdir, 'dir')
    mkdir(outdir);
end

scl = 0.25;
spacing = 75 * scl;
se = strel('disk', round(spacing));

ext = '.png';

imglist = dir([imgpath, filesep, '*', ext]);
discardlist = dir([discardpath, filesep, '*', ext]);
discardidx = ismember({imglist.name}, {discardlist.name});
imglist = imglist(~discardidx);
allnames = {imglist.name};
Nimg = numel(imglist);

covered_pt = cell(Nimg, 2);

matlist = dir([ptpairpath, filesep, '*.mat']);

for k = 1:numel(matlist)
    load([matlist(k).folder, filesep, matlist(k).name], 'imgnames');
    [~, sec_ids] = ismember(imgnames, allnames);
    if any(sec_ids == 0)
        continue
    end
    idx1 = sec_ids(1);
    idx2 = sec_ids(2);
    load([matlist(k).folder, filesep, matlist(k).name], 'ptpairs');
    yx1 = ptpairs.yx1;
    yx2 = ptpairs.yx2;
    if idx1 < idx2
        covered_pt{idx1,1} = [covered_pt{idx1,1}; yx1];
        covered_pt{idx2,2} = [covered_pt{idx2,2}; yx2];
    else
        covered_pt{idx1,2} = [covered_pt{idx1,2}; yx1];
        covered_pt{idx2,1} = [covered_pt{idx2,1}; yx2];
    end
end
%%
for k = 1:Nimg
    if contains(imglist(k).name, 'stitch')
        continue
    end
    IMG = imread([imglist(k).folder, filesep, imglist(k).name]);
    IMGt = single(imresize(IMG, scl));
    imgsz = size(IMGt);
    yx1 = covered_pt{k, 1};
    yx2 = covered_pt{k, 2};
    mask1 = zeros(imgsz, 'single');
    mask2 = zeros(imgsz, 'single');
    if ~isempty(yx1)
        yx1 = min(max(1, round(yx1 * scl)), imgsz);
        indx = sub2ind(imgsz, yx1(:,1), yx1(:,2));
        mask1(indx) = 1;
        mask1 = imdilate(mask1, se);
    end
    if ~isempty(yx2)
        yx2 = min(max(1, round(yx2 * scl)), imgsz);
        indx = sub2ind(imgsz, yx2(:,1), yx2(:,2));
        mask2(indx) = 1;
        mask2 = imdilate(mask2, se);
    end
    IMGo = uint8(cat(3, IMGt, IMGt .* mask1, IMGt .* mask2));
    imwrite(IMGo, [outdir, filesep, imglist(k).name]);
end
