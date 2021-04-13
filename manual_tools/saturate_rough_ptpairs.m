clear; clc;
addpath('workflows','align_config')
parentpath = get_parent_path;

gridsz = 300;
ext = '.png';


imgpath = [parentpath, filesep, 'stitched'];
maskrpath = [parentpath, filesep, 'masks_region'];
maskwpath = [parentpath, filesep, 'masks_wrinkle'];
outpath = [parentpath, filesep, 'ptpairs0'];
imglist = dir([imgpath, filesep, '*', ext]);
matlist = dir([outpath, filesep, '*.mat']);

pt_idxs = 1:numel(matlist);% find(strcmpi({matlist.name},'0485_0486.mat'));

Nimg = numel(imglist);
if isempty(imglist)
    return
end

if ~exist(outpath, 'dir')
    mkdir(outpath);
end
warning('off')
t0 = 0;
tic;

for p = 1:numel(pt_idxs)
    kp = pt_idxs(p);
    outname = [matlist(kp).folder, filesep, matlist(kp).name];
    try
        load(outname);
        disp(outname)
        A = geometries.fit_affine(ptpairs.yx2, ptpairs.yx1);
        [U,S,V] = svd(A(1:2, 1:2));
        R = U * V';
        rot = atan2(R(2), R(1));
        maskr0 = imread([maskrpath, filesep, imgnames{1}]);
        maskw0 = imread([maskwpath, filesep, imgnames{1}]);
        maskr1 = imread([maskrpath, filesep, imgnames{2}]);
        maskw1 = imread([maskwpath, filesep, imgnames{2}]);

        maskw0 = imdilate(imdilate(maskw0 == 0, ones(gridsz,1)), ones(1,gridsz));
        maskw1 = imdilate(imdilate(maskw1 == 0, ones(gridsz,1)), ones(1,gridsz));

        [gx, gy] = meshgrid(1:gridsz:size(maskr0,2), 1:gridsz:size(maskr0,1));
        gyx0 = round([gy(:), gx(:)]);
        gyx1 = round(gyx0 * A(1:2,1:2) + A(3,1:2));
        valid_idx = all((gyx0 > 0) & (gyx0 < size(maskr0)) & (gyx1 > 0) & (gyx1 < size(maskr1)),2);
        gyx0 = gyx0(valid_idx,:);
        gyx1 = gyx1(valid_idx,:);
        gindx0 = sub2ind(size(maskr0), gyx0(:,1), gyx0(:,2));
        gindx1 = sub2ind(size(maskr1), gyx1(:,1), gyx1(:,2));
        rid0 = maskr0(gindx0);
        rid1 = maskr1(gindx1);
        wid0 = maskw0(gindx0);
        wid1 = maskw1(gindx1);
        valid_idx = (rid0>0) & (rid1>0) & wid0 & wid1;
        yx1 = gyx0(valid_idx, :);
        yx2 = gyx1(valid_idx, :);
        region_id = [rid0(valid_idx), rid1(valid_idx)];
        conf = ones(sum(valid_idx),1, 'single');
        rotation = rot * ones(sum(valid_idx),1, 'single');
        ptpairs = table(yx1,yx2,region_id,conf, rotation);
        save(outname, 'ptpairs', 'imgnames');
        shr = exp(abs(log(diag(S))));
        if  max(shr) > 1.25
            disp([matlist(kp).name, ': shear ', num2str(max(shr))]);
        end
    catch ME
        disp(matlist(kp).name);
        disp(ME.message);
    end
    t0 = t0 + 1;
    if mod(t0, 30) == 0
        tt = toc;
        disp([num2str(tt / t0),'sec']);
    end
end