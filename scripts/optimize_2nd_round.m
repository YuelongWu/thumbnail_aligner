clear; clc;
addpath('workflows', 'align_config')
use_parallel = false;
startover = false;
parentpath = get_parent_path;
options = parse_yaml(['align_config', filesep, 'optimize.yaml']);
ext = '.png';
imgpath = [parentpath, filesep, 'stitched'];
secondpath = [imgpath, filesep, '2nd_round'];
searchpaths = {imgpath, secondpath};

preylist = dir([secondpath, filesep, '*.png']);
% preylist = dir([imgpath, filesep,'*.png']);
preyname = erase({preylist.name}, ext);
% preyname = preyname(contains(preyname, 'stitch3'));
% preyname = {'sec2240ref_stitch3'};

ptpairpath = [parentpath, filesep, 'ptpairs'];
maskrpath = [parentpath, filesep, 'masks_region'];
maskwpath = [parentpath, filesep, 'masks_wrinkle'];
meshpath = [parentpath, filesep, 'meshes'];

mlist = dir([meshpath, filesep, '*.mat']);
if ~isempty(mlist)
    processed_name = erase({mlist.name}, '.mat');
    all_relevant_name = unique([processed_name, preyname]);
    load([mlist(1).folder, filesep, mlist(1).name], 'outsz');
    if ~startover
        idxt = ismember(preyname, processed_name);
        preyname = preyname(~idxt);
    end
else
    outimgpath = [parentpath, filesep, 'output'];
    mlist = dir([outimgpath, filesep, '*', ext]);
    processed_name = erase({mlist.name}, ext);
    all_relevant_name = unique([processed_name, preyname]);
    IMG0 = imread([mlist(1).folder, filesep, mlist(1).name]);
    outsz = size(IMG0);
end

matlist = dir([ptpairpath, filesep, '*.mat']);
[idx, isprey] = utils.filter_img_from_mat(all_relevant_name, {matlist.name}, preyname);
all_relevant_name = all_relevant_name(idx);
fixed_sec = ~isprey(:);

N_all = numel(all_relevant_name);
imglist0 = cell(N_all, 1);
for k = 1:N_all
    [~, imgdir] = utils.try_read_img(searchpaths, [all_relevant_name{k}, '.png'], false);
    if isempty(imgdir)
        error(['missing image: ', all_relevant_name{k}]);
    end
    imglist0{k} = imgdir;
end
imglist0 = vertcat(imglist0{:});

%%
Nprey = numel(preyname);
imglistC = cell(Nprey,1);
fixed_secC = cell(Nprey,1);
pplistC = cell(Nprey,1);
MATnames = erase({matlist.name},'.mat');
MATnames_split = utils.split_pair_name(MATnames(:));
IMGnames = erase({imglist0.name}, ext);
IMGnames = IMGnames(:);
prey_allocated = false(Nprey, 1);
for k = 1:numel(preyname)
    if prey_allocated(k)
        continue
    end
    idxtPREY = strcmpi(preyname, preyname{k});
    idxtIMG = strcmpi(IMGnames, preyname{k});
    cnt1 = sum(idxtPREY);
    cnt0 = 0;
    while(sum(idxtIMG) > cnt0)
        cnt0 = sum(idxtIMG);
        idxt = utils.filter_img_from_mat(IMGnames, MATnames, preyname(idxtPREY));
        idxtIMG(idxt) = true;
        idxtPREY = ismember(preyname, IMGnames(idxtIMG));
        prey_allocated(idxtPREY) = true;
        if sum(idxtPREY) == cnt1
            break
        else
            cnt1 = sum(idxtPREY);
        end
    end
    imglistC{k} = imglist0(idxtIMG);
    fixed_secC{k} = fixed_sec(idxtIMG);
    idxtPP = ismember(MATnames_split(:,1), preyname(idxtPREY)) | ismember(MATnames_split(:,2), preyname(idxtPREY));
    pplistC{k} = matlist(idxtPP);
end
idxt = cellfun(@isempty, imglistC);
imglistC = imglistC(~idxt);
pplistC = pplistC(~idxt);
fixed_secC = fixed_secC(~idxt);

Npool = 5;
Nthrd = 4;
if use_parallel
    pool_obj = gcp('nocreate');
    if isempty(pool_obj) || ~pool_obj.Connected
        pool_obj = parpool(Npool);
    else
        Npool = pool_obj.NumWorkers;
    end
end

%%
optop = options.optimizer;
for f0 = 1:numel(imglistC)
    tt0 = tic;
    imglist_f = imglistC{f0};
    matlist_f = pplistC{f0};
    fixed_sec_f = fixed_secC{f0};
%     allexist = utils.check_all_file_exists(meshpath, strrep({imglist_f.name},ext,'.mat'));
%     if allexist
%         continue
%     end
    [Ms, cost_history] = optimize_wrapper_multires(imglist_f, matlist_f, fixed_sec_f, ...
        optop, maskwpath, maskrpath, meshpath);

    for k = 1:numel(Ms)
        outname = [meshpath, filesep, strrep(imglist_f(k).name, ext, '.mat')];
        if fixed_sec_f(k) && exist(outname, 'file')
            continue;
        end
        % if exist([meshpath, filesep, outname],'file')
        %     continue;
        % end
        M = Ms{k};
        if exist([maskrpath, filesep, imglist_f(k).name], 'file')
            maskr = imread([maskrpath, filesep, imglist_f(k).name]);
        else
            img = imread([imglist_f(k).folder, filesep, imglist_f(k).name]);
            maskr = ones(size(img), 'uint8');
        end
        if exist([maskwpath, filesep, imglist_f(k).name], 'file')
            maskw = imread([maskwpath, filesep, imglist_f(k).name]);
        else
            img = imread([imglist_f(k).folder, filesep, imglist_f(k).name]);
            maskw = zeros(size(img), 'uint8');
        end
        if use_parallel
            utils.parsave(outname, {M, maskr, maskw, outsz}, {'M', 'maskr','maskw', 'outsz'});
        else
            save(outname, 'M', 'maskr','maskw', 'outsz');
        end
    end
    toc(tt0);
end

%%
% pause(300);
% block_tform_render_whole_section