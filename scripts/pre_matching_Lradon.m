clear; clc;
addpath('workflows','align_config')
parentpath = get_parent_path;
options = parse_yaml(['align_config', filesep, 'subregion_matching.yaml']);
Nstep = 3;
ext = '.png';
imgpath = [parentpath, filesep, 'stitched'];
maskrpath = [parentpath, filesep, 'masks_region'];
maskwpath = [parentpath, filesep, 'masks_wrinkle'];
outpath = [parentpath, filesep, 'ptpairs0'];

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
            % if abs(str2double(erase(mname0,{'ref', 'sec'})) - str2double(erase(mname1,{'ref', 'sec'}))) > 3
            %     continue
            % end
            outname = [mname0,'_', mname1,'.mat'];
            % outname = [strrep(imglist(k0).name, ext, ''),'_', strrep(imglist(k1).name, ext, ''),'.mat'];
            if exist([outpath, filesep, outname], 'file')
                continue
            end
            % if any(kk ~= [123, 124])
            %     continue
            % end
            if any(des_ids == k0)
                info0 = des_stack{find(des_ids == k0,1)};
            else
                IMG0 = imread([imglist(k0).folder, filesep, imglist(k0).name]);
                if exist([maskrpath, filesep, imglist(k0).name], 'file')
                    maskr0 = imread([maskrpath, filesep, imglist(k0).name]);
                else
                    maskr0 = ones(size(IMG0), 'uint8');
                end
                if exist([maskwpath, filesep, imglist(k0).name], 'file')
                    maskw0 = imread([maskwpath, filesep, imglist(k0).name]);
                else
                    maskw0 = zeros(size(IMG0), 'uint8');
                end
                info0 =  local_radon.raw_img_to_descriptor_subregion(IMG0, options, {maskw0, maskr0});
                idxt = find(isnan(des_ids), 1);
                des_stack{idxt} = info0;
                des_ids(idxt) = k0;
            end
            if any(des_ids == k1)
                info1 = des_stack{find(des_ids == k1,1)};
            else
                IMG1 = imread([imglist(k1).folder, filesep, imglist(k1).name]);
                if exist([maskrpath, filesep, imglist(k1).name], 'file')
                    maskr1 = imread([maskrpath, filesep, imglist(k1).name]);
                else
                    maskr1 = ones(size(IMG1), 'uint8');
                end
                if exist([maskwpath, filesep, imglist(k1).name], 'file')
                    maskw1 = imread([maskwpath, filesep, imglist(k1).name]);
                else
                    maskw1 = zeros(size(IMG1), 'uint8');
                end
                info1 =  local_radon.raw_img_to_descriptor_subregion(IMG1, options, {maskw1, maskr1});
                idxt = find(isnan(des_ids), 1);
                des_stack{idxt} = info1;
                des_ids(idxt) = k1;
            end
            [ptpairs, cnt, ratio, svs, rgpairs] = subregion_matching(info0, info1, {maskw0, maskw1}, {maskr0, maskr1}, options);
            imgnames = {imglist(k0).name, imglist(k1).name};
            save([outpath, filesep, outname], 'ptpairs', 'cnt', 'ratio', 'svs', 'rgpairs', 'kk','imgnames');
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