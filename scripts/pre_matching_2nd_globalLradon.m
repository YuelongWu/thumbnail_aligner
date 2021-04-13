clear; clc;
addpath('workflows','align_config')
use_parallel = false;
parentpath = get_parent_path;
options = parse_yaml(['align_config', filesep, 'initial_matching.yaml']);
Nstep = [-1,1];
scl = 1;
ext = '.png';
imgpath = [parentpath, filesep, 'stitched'];
maskrpath = [parentpath, filesep, 'masks_region'];
maskwpath = [parentpath, filesep, 'masks_wrinkle'];
outpath = [parentpath, filesep, 'ptpairs0'];
secondpath = [imgpath, filesep, '2nd_round'];

imglist = dir([imgpath, filesep, '*', ext]);
secondlist = dir([secondpath, filesep, '*', ext]);

issecond = [false(numel(imglist),1); true(numel(secondlist), 1)];
imglist = [imglist(:);secondlist(:)];
if isempty(imglist) || isempty(secondlist)
    return
end
IMGnames = {imglist.name};
[~, idxt] = sort(IMGnames);
imglist = imglist(idxt);
issecond = issecond(idxt);

Nsec = numel(secondlist);


idx_second = find(issecond);
secondlist = imglist(idx_second);
secdis_fwd = cumsum(1-issecond(:));
secdis_bwd = -flipud(cumsum(flipud(1-issecond(:))));

partnerids = cell(Nsec, 1);
partnerlists = cell(Nsec, 1);
for k = 1:Nsec
    idxt = idx_second(k);
    disfwd_t = max(0, secdis_fwd - secdis_fwd(idxt));
    disbwd_t = min(0, secdis_bwd - secdis_bwd(idxt));
    secdis_t = (disfwd_t + disbwd_t) .* (~issecond);
    [~, ib] = ismember(Nstep(:), secdis_t(:));
    partnerids{k} = ib(ib > 0);
    partnerlists{k} = imglist(ib(ib > 0));
end

Nimg = numel(imglist);

if ~exist(outpath, 'dir')
    mkdir(outpath);
end

Npool = 5;
Nthrd = 2;
if use_parallel
    pool_obj = gcp('nocreate');
    if isempty(pool_obj) || ~pool_obj.Connected
        pool_obj = parpool(Npool);
    else
        Npool = pool_obj.NumWorkers;
    end
end

warning('off');
for f0 = 1:Nsec
    tic
    info0_created = false;
    k0 = idx_second(f0);
    ptid = partnerids{f0};
    if isempty(ptid)
        continue
    end
    ptlist = partnerlists{f0};
    for f1 = 1:numel(ptid)
        k1 = ptid(f1);
        kk = [k0, k1];
        try
            outname = [strrep(secondlist(f0).name, ext, ''),'_', strrep(ptlist(f1).name, ext, ''),'.mat'];
            if exist([outpath, filesep, outname], 'file')
                continue
            end
            if ~info0_created
                IMG0 = imread([secondlist(f0).folder, filesep, secondlist(f0).name]);
                if exist([maskrpath, filesep, secondlist(f0).name], 'file')
                    maskr0 = imread([maskrpath, filesep, secondlist(f0).name]);
                    maskr0 = uint8(maskr0 > 0);
                else
                    maskr0 = ones(size(IMG0), 'uint8');
                end
                if exist([maskwpath, filesep, secondlist(f0).name], 'file')
                    maskw0 = imread([maskwpath, filesep, secondlist(f0).name]);
                else
                    maskw0 = zeros(size(IMG0), 'uint8');
                end
                if scl < 1
                    IMG0 = imresize(IMG0, scl, 'bilinear');
                    maskr0 = imresize(maskr0, scl, 'nearest');
                    maskw0 = imresize(maskw0, scl, 'nearest');
                end
                info0 =  local_radon.raw_img_to_descriptor_subregion(IMG0, options, {maskw0, maskr0});
                info0_created = true;
            end

            IMG1 = imread([ptlist(f1).folder, filesep, ptlist(f1).name]);
            if exist([maskrpath, filesep, ptlist(f1).name], 'file')
                maskr1 = imread([maskrpath, filesep, ptlist(f1).name]);
                maskr1 = uint8(maskr1 > 0);
            else
                maskr1 = ones(size(IMG1), 'uint8');
            end
            if exist([maskwpath, filesep, ptlist(f1).name], 'file')
                maskw1 = imread([maskwpath, filesep, ptlist(f1).name]);
            else
                maskw1 = zeros(size(IMG1),'uint8');
            end
            if scl < 1
                IMG1 = imresize(IMG1, scl, 'bilinear');
                maskr1 = imresize(maskr1, scl, 'nearest');
                maskw1 = imresize(maskw1, scl, 'nearest');
            end
            info1 =  local_radon.raw_img_to_descriptor_subregion(IMG1, options, {maskw1, maskr1});

            [ptpairs, Np] = global_matching(info0, info1,options);
            imgnames = {secondlist(f0).name, ptlist(f1).name};
            if scl < 1
                ptpairs.yx1 = ptpairs.yx1 / scl;
                ptpairs.yx2 = ptpairs.yx2 / scl;
            end
            if use_parallel
                utils.parsave([outpath, filesep, outname], {ptpairs, imgnames}, {'ptpairs','imgnames'});
            else
                save([outpath, filesep, outname], 'ptpairs', 'imgnames');
            end

        catch ME
            disp([k0, k1]);
            disp(ME.message);
        end
    end
    toc;
end
