clear; clc;
addpath('workflows','align_config')
use_parallel = true;
parentpath = get_parent_path;
options = parse_yaml(['align_config', filesep, 'initial_matching.yaml']);

lstep = [1,2,5];
Nstep = numel(lstep);
scl = 0.5;

ext = '.png';
imgpath = [parentpath, filesep, 'stitched'];
maskrpath = [parentpath, filesep, 'masks_region'];
maskwpath = [parentpath, filesep, 'masks_wrinkle'];
outpath = [parentpath, filesep, 'ptpairs0'];

if exist([parentpath, filesep, 'imglist.txt'],'file')
    fid = fopen([parentpath, filesep, 'imglist.txt']);
    outcell = textscan(fid, '%s');
    fclose(fid)
    outcell = outcell{1};
    imglist = struct('folder', imgpath, 'name', strcat(outcell,'.png'));
else
    imglist = dir([imgpath, filesep, '*', ext]);
end
Nimg = numel(imglist);
if isempty(imglist)
    return
end

if ~exist(outpath, 'dir')
    mkdir(outpath);
end

Npool = 3;
Nthrd = max(1,floor(6/Npool));
if use_parallel
    pool_obj = gcp('nocreate');
    if isempty(pool_obj) || ~pool_obj.Connected
        pool_obj = parpool(Npool);
    else
        Npool = pool_obj.NumWorkers;
    end
end
Njob = min(Npool * 3, Nimg);
idxjob = round(linspace(1, Nimg, Njob+1));

warning('off');
parfor kjob = 1:Njob
idx0 = idxjob(kjob);
idx1 = idxjob(kjob+1);
des_stack = cell(Nstep + 2, 1);
des_ids = nan(Nstep + 2, 1);
des_cnt = inf(Nstep + 2, 1);
t0 = 0;
tic;
for k0 = idx0:idx1
    for k1t = 1:Nstep
        k1 = k0 + lstep(k1t);
        if k1 > Nimg
            continue
        end
        kk = [k0, k1];
        try
            mname0 = erase(imglist(k0).name, ext);
            mname1 = erase(imglist(k1).name, ext);
            if contains(mname0, 'W19') && contains(mname1, 'W19')
                continue
            end
            % if contains(mname0, 'stitch3')
            %     continue
            % end
            % if ~strcmpi(erase(mname0, '_goog14'), erase(mname1, {'_stitch3','R1','R2','R3'}))
            %     continue
            % end
            % if abs(str2double(erase(mname0,{'ref', 'sec'})) - str2double(erase(mname1,{'ref', 'sec'}))) > 3
            %     continue
            % end
            % if k0 >= k1
            %     continue
            % end
            outname = [mname0,'_', mname1,'.mat'];
            % save([outpath, filesep, outname],'outname');
            if exist([outpath, filesep, outname], 'file')
                continue
            end
            % if any(kk ~= [123, 124])
            %     continue
            % end
            if any(des_ids == k0)
                info0 = des_stack{find(des_ids == k0,1)};
                des_cnt(des_ids==k0) = 1;
            else
                IMG0 = imread([imglist(k0).folder, filesep, imglist(k0).name]);
                if exist([maskrpath, filesep, imglist(k0).name], 'file')
                    maskr0 = imread([maskrpath, filesep, imglist(k0).name]);
                else
                    maskr0 = ones(size(IMG0), 'uint8');
                end
                maskr0 = uint8(maskr0 > 0);
                if exist([maskwpath, filesep, imglist(k0).name], 'file')
                    maskw0 = imread([maskwpath, filesep, imglist(k0).name]);
                else
                    maskw0 = zeros(size(IMG0), 'uint8');
                end
                if scl < 1
                    IMG0 = imresize(IMG0, scl, 'bilinear');
                    maskr0 = imresize(maskr0, scl, 'nearest');
                    maskw0 = imresize(maskw0, scl, 'nearest');
                end
                info0 =  local_radon.raw_img_to_descriptor_subregion(IMG0, options, {maskw0, maskr0});
                [~, idxt] = max(des_cnt);
                des_cnt = des_cnt + 1;
                des_stack{idxt} = info0;
                des_ids(idxt) = k0;
                des_cnt(idxt) = 1;
            end
            if any(des_ids == k1)
                info1 = des_stack{find(des_ids == k1,1)};
                des_cnt(des_ids==k1) = 1;
            else
                IMG1 = imread([imglist(k1).folder, filesep, imglist(k1).name]);
                if exist([maskrpath, filesep, imglist(k1).name],'file')
                    maskr1 = imread([maskrpath, filesep, imglist(k1).name]);
                else
                    maskr1 = ones(size(IMG1),'uint8');
                end
                maskr1 = uint8(maskr1 > 0);
                if exist([maskwpath, filesep, imglist(k1).name],'file')
                    maskw1 = imread([maskwpath, filesep, imglist(k1).name]);
                else
                    maskw1 = zeros(size(IMG1),'uint8');
                end
                if scl < 1
                    IMG1 = imresize(IMG1, scl, 'bilinear');
                    maskr1 = imresize(maskr1, scl, 'nearest');
                    maskw1 = imresize(maskw1, scl, 'nearest');
                end
                info1 =  local_radon.raw_img_to_descriptor_subregion(IMG1, options, {maskw1, maskr1});
                [~,idxt] = max(des_cnt);
                des_cnt = des_cnt + 1;
                des_stack{idxt} = info1;
                des_ids(idxt) = k1;
                des_cnt(idxt) = 1;
            end
            [ptpairs, Np] = global_matching(info0, info1, options);
            imgnames = {imglist(k0).name, imglist(k1).name};
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
            disp(outname);
            disp(ME.message);
        end
        t0 = t0 + 1;
        if mod(t0, 30) == 0
            tt = toc;
            disp([num2str(tt / t0),'sec']);
        end
    end
    % des_ids(des_ids == k0) = nan;
end
end