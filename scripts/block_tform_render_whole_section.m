clc;clear
addpath('align_config');
options = parse_yaml(['align_config', filesep, 'optimize.yaml']);
render_option = options.render;
parentpath = get_parent_path;
ext = '.png';

scl = 1;
intp = 'linear';
outputall = true;
imgpath = [parentpath, filesep, 'stitched'];
secondpath = [imgpath, filesep, '2nd_round'];
imgpaths = {imgpath, secondpath};
fillval = 255;

use_parallel = true;
Npool = 2;
Nthrd = 3;

meshpath = [parentpath, filesep, 'meshes'];
outimgpath = [parentpath, filesep, 'output'];
outxypath = [parentpath, filesep, 'displacement'];

if ~exist(outimgpath, 'dir')
    mkdir(outimgpath);
end
if ~exist(outxypath, 'dir')
    mkdir(outxypath);
end

matlist = dir([meshpath, filesep, '*.mat']);
if use_parallel
    pool_obj = gcp('nocreate');
    if isempty(pool_obj) || ~pool_obj.Connected
        pool_obj = parpool(Npool);
    else
        Npool = pool_obj.NumWorkers;
    end
end

parfor k = 1:numel(matlist)                                                     % par
    if use_parallel
        maxNumCompThreads(Nthrd);                                                % par
    end
    % maxNumCompThreads(Nthrd);  
    try
        imgname = strrep(matlist(k).name, '.mat', ext);
        % if contains(imgname, 'goog14')
        %     continue
        % end
        if exist([outimgpath, filesep, imgname], 'file')
            continue
        end
        tic
        meshdir = [matlist(k).folder, filesep, matlist(k).name];
        M = utils.parload(meshdir, 'M');
        maskr = utils.parload(meshdir, 'maskr');
        maskw = utils.parload(meshdir, 'maskw');
        maskw = uint8(maskw == 200);
        outsz = utils.parload(meshdir, 'outsz');
        masks = cat(3, maskr, maskw);
        MC = elastic_mesh.mesh2mov_block(M, masks, render_option, [], outsz, outputall);
        IMG = utils.try_read_img(imgpaths, imgname);
        [IMGt, D]= elastic_mesh.block_tform_render_whole(IMG, MC, scl, intp, fillval);
        if any(size(IMGt) ~= outsz)
            IMGt0 = fillval * ones(outsz, class(IMGt));
            imght = min(size(IMGt0,1), size(IMGt,1));
            imgwd = min(size(IMGt0,2), size(IMGt,2));
            IMGt0(1:imght, 1:imgwd) = IMGt(1:imght, 1:imgwd);
            IMGt = IMGt0;
            D0 = zeros(outsz(1),outsz(2),2,class(D));
            imght = min(size(D0,1), size(D,1));
            imgwd = min(size(D0,2), size(D,2));
            D0(1:imght, 1:imgwd,:) = D(1:imght, 1:imgwd,:);
            D = D0;
        end
        imwrite(IMGt, [outimgpath, filesep, imgname]);
        if M.region_num > 1
            utils.parsave([outxypath, filesep, matlist(k).name], D, 'D');
        end
        t = toc;
        disp(['section ', imgname,'| worker ', num2str(feature('getpid')), '| ', num2str(t, '%5f'), ' sec']);
    catch ME
        disp(['error: ', imgname])
        disp(ME)
    end
end

% delete(poolobj);                                                           % par