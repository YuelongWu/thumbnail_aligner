clear; clc;

addpath('workflows','align_config')
parentpath = get_parent_path;
% options = parse_yaml(['align_config', filesep, 'subregion_matching.yaml']);
Nstep = 4;
gridsz = 200;
ext = '.png';
imgpath = [parentpath, filesep, 'stitched'];
maskrpath = [parentpath, filesep, 'masks_region'];
maskwpath = [parentpath, filesep, 'masks_wrinkle'];
outpath = [parentpath, filesep, 'ptpairs0_dummpy'];

imglist = dir([imgpath, filesep, '*', ext]);
Nimg = numel(imglist);
if isempty(imglist)
    return
end

if ~exist(outpath, 'dir')
    mkdir(outpath);
end
warning('off')
des_stack = cell(Nstep + 2, 1);
des_ids = nan(Nstep + 2, 1);
t0 = 0;
tic;

for k0 = 1:Nimg
    for k1 = (k0 + 1):min(k0 + Nstep, Nimg)
        kk = [k0, k1];
        try
            mname0 = erase(imglist(k0).name, ext);
            mname1 = erase(imglist(k1).name, ext);
            if contains(mname0,'stitch3')
                continue
            end
            if abs(str2double(erase(mname0,{'ref', 'sec'})) - str2double(erase(mname1,{'ref', 'sec','stitch3'}))) > 0
                continue
            end
            outname = [mname0,'_', mname1,'.mat'];
            if exist([outpath, filesep, outname], 'file')
               continue
            end
            % if any(kk ~= [123, 124])
            %     continue
            % end
            maskr0 = imread([maskrpath, filesep, imglist(k0).name]);
            % maskw0 = imread([maskwpath, filesep, imglist(k0).name]);
            maskr1 = imread([maskrpath, filesep, imglist(k1).name]);
            % maskw1 = imread([maskwpath, filesep, imglist(k1).name]);
            [gx, gy] = meshgrid(1:gridsz:size(maskr0,2), 1:gridsz:size(maskr0,1));
            gyx0 = [gy(:), gx(:)];
            gyx1 = gyx0;
            valid_idx = all((gyx0 > 0) & (gyx0 < size(maskr0)) & (gyx1 > 0) & (gyx1 < size(maskr1)),2);
            gyx0 = gyx0(valid_idx,:);
            gyx1 = gyx1(valid_idx,:);
            gindx0 = sub2ind(size(maskr0), gyx0(:,1), gyx0(:,2));
            gindx1 = sub2ind(size(maskr1), gyx1(:,1), gyx1(:,2));
            rid0 = maskr0(gindx0);
            rid1 = maskr1(gindx1);
            valid_idx = (rid0>0) & (rid1>0);
            yx1 = gyx0(valid_idx, :);
            yx2 = gyx1(valid_idx, :);
            region_id = [rid0(valid_idx), rid1(valid_idx)];
            conf = ones(sum(valid_idx),1);
            rotation = zeros(sum(valid_idx),1);
            ptpairs = table(yx1,yx2,region_id,conf, rotation);
            [rgpairs,~,ic] = unique(region_id, 'rows');
            cnt = histcounts(ic, 0.5:1:(max(ic(:))+0.5)); cnt = cnt(:);
            ratio = ones(size(cnt));
            svs = zeros(size(cnt));
            imgnames = {imglist(k0).name, imglist(k1).name};
            save([outpath, filesep, outname], 'ptpairs','imgnames');
            if max(svs(ratio > 0.25)) > 0.1 || min(cnt(ratio > 0.25) < 5)
                disp(['missing_match:',imgnames{1}, '->', imgnames{2}])
            end
        catch ME
            disp([k0, k1]);
            disp(ME.message);
        end
        t0 = t0 + 1;
        if mod(t0, 30) == 0
            tt = toc;
            disp([num2str(tt / t0),'sec']);
        end
    end
    des_ids(des_ids == k0) = nan;
end