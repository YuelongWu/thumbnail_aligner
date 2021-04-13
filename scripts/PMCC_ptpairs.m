clc;clear;
use_parallel = true;
addpath('workflows','align_config');
parentpath = get_parent_path;
% options =  parse_yaml('E:\one_time_stuff\jboulanger\google_pipeline\prematch.yaml');
% options = options.finematch;
options = parse_yaml(['align_config', filesep, 'PMCC_iter.yaml']);
kfpath = [parentpath, filesep, 'keyFrames'];
ptpairpath = [parentpath, filesep, 'ptpairs0'];
ext = '.png';
imgpath = [parentpath, filesep, 'stitched'];
secondpath = [imgpath, filesep, '2nd_round'];
workpath = [imgpath, filesep, 'recently_added'];
maskrpath = [parentpath, filesep, 'masks_region'];
maskwpath = [parentpath, filesep, 'masks_wrinkle'];
% outpath = [kfpath, filesep, 'ptpairs'];
imgpaths = {imgpath, secondpath, workpath};

matlist0 = dir([ptpairpath, filesep, '*.mat']);
% matlist0 = dir([kfpath,filesep,'**',filesep,'ptpairs0', filesep,'*.mat']);
if isempty(matlist0)
    return
end
Nmat = numel(matlist0);

Npool = 3;
Nthrd = 2;
if use_parallel
    pool_obj = gcp('nocreate');
    if isempty(pool_obj) || ~pool_obj.Connected
        pool_obj = parpool(Npool);
    else
        Npool = pool_obj.NumWorkers;
    end
end
Njob = min(Npool * 5, Nmat);
matlistCell = cell(Njob, 1);
idxt = round(linspace(1, Nmat, numel(matlistCell)+1));
for k = 1:numel(matlistCell)
    matlistCell{k} = matlist0(idxt(k):idxt(k+1));
end

% if ~exist(outpath, 'dir')
%     mkdir(outpath);
% end
warning('off');

parfor kbtch = 1:Njob
matlist = matlistCell{kbtch};
if use_parallel
    maxNumCompThreads(Nthrd);
else
    maxNumCompThreads(Nthrd * Npool);
end
cache_sz = 5;
cache_info = cell(cache_sz,1);
cache_names = cell(cache_sz,1);
t0 = 0;
ctt = 1;
tic;
for k = 1:1:numel(matlist)
    try
        % if ~strcmpi(matlist(k).name,'sec2641_sec2642ref.mat')
        %     continue
        % end
        outpath = strrep([matlist(k).folder], 'ptpairs0','ptpairs');
        if ~exist(outpath, 'dir')
            mkdir(outpath);
        end
        if exist([outpath, filesep, matlist(k).name], 'file')
            continue
        end
        if use_parallel
            ptpairs = utils.parload([matlist(k).folder, filesep, matlist(k).name], 'ptpairs');
            imgnames = utils.parload([matlist(k).folder, filesep, matlist(k).name], 'imgnames');
        else
            load([matlist(k).folder, filesep, matlist(k).name], 'ptpairs', 'imgnames');
        end
        ptpairs0 = ptpairs;
        if any(strcmpi(cache_names, imgnames{1}))
            info0 = cache_info{find(strcmpi(cache_names, imgnames{1}),1)};
            maskr0 = info0.maskr;
            maskw0 = info0.maskw;
        else
            % IMG0 = imread([imgpath, filesep, imgnames{1}]);
            IMG0 = utils.try_read_img(imgpaths, imgnames{1});
            try
                maskr0 = imread([maskrpath, filesep, imgnames{1}]);
            catch
                maskr0 = ones(size(IMG0), 'uint8');
            end
            try
                maskw0 = imread([maskwpath, filesep, imgnames{1}]);
            catch
                maskw0 = zeros(size(IMG0), 'uint8');
            end
            % info0 =  PMCC.raw_img_to_descriptor(IMG0, options, {maskw0, maskr0});
            % cache_info{mod(t0, cache_sz) + 1} = info0;
            % cache_names{mod(t0, cache_sz) + 1} = imgnames{1};
            % t0 = t0 + 1;
            info0 = IMG0;
        end

        if any(strcmpi(cache_names, imgnames{2}))
            info1 = cache_info{find(strcmpi(cache_names, imgnames{2}),1)};
            maskr1 = info1.maskr;
            maskw1 = info1.maskw;
        else
            % IMG1 = imread([imgpath, filesep, imgnames{2}]);
            IMG1 = utils.try_read_img(imgpaths, imgnames{2});
            try
                maskr1 = imread([maskrpath, filesep, imgnames{2}]);
                maskr1([1,end],:) = 0;
                maskr1(:,[1,end]) = 0;
            catch
                maskr1 = ones(size(IMG1), 'uint8');
            end
            try
                maskw1 = imread([maskwpath, filesep, imgnames{2}]);
            catch
                maskw1 = zeros(size(IMG1), 'uint8');
            end
            % info1 =  PMCC.raw_img_to_descriptor(IMG1, options, {maskw1, maskr1});
            % cache_info{mod(t0, cache_sz) + 1} = info1;
            % cache_names{mod(t0, cache_sz) + 1} = imgnames{2};
            % t0 = t0 + 1;
            info1 = IMG1;
        end
        [ptpairs, info0t, info1t] = PMCC_matching_iter(info0, info1, {maskw0, maskw1}, {maskr0, maskr1}, options, ptpairs0);
        if isstruct(info0t) && ~isstruct(info0)
            cache_info{mod(t0, cache_sz) + 1} = info0t;
            cache_names{mod(t0, cache_sz) + 1} = imgnames{1};
            t0 = t0 + 1;
        end
        if isstruct(info1t) && ~isstruct(info1)
            cache_info{mod(t0, cache_sz) + 1} = info1t;
            cache_names{mod(t0, cache_sz) + 1} = imgnames{2};
            t0 = t0 + 1;
        end
        % ptpairs = [ptpairs0;ptpairs];
        ptpairs.conf = ones(size(ptpairs.conf));
        % ptpairs.conf = max(0, ptpairs.conf/max(ptpairs.conf(:))).^0.5;
        if use_parallel
            utils.parsave([outpath, filesep, matlist(k).name], {ptpairs, imgnames}, {'ptpairs','imgnames'});
        else
            save([outpath, filesep, matlist(k).name],'ptpairs','imgnames');
        end
        ctt = ctt + 1;
    catch ME
        disp(matlist(k).name);
        disp(ME.message)
    end
    if mod(ctt, 10) == 0
        tt = toc;
        disp([num2str(tt/ctt), 'sec']);
        tic;
        ctt = 0;
    end
end
end
