clc;clear;
addpath('align_config', 'workflows');
options = parse_yaml(['align_config', filesep, 'optimize.yaml']);
parentpath = get_parent_path;
ext = '.png';

shareM0 = false;
rigidFirst = true;
affineFirst = true;
toRender = false;
elastic = true;
full_cover = true;
optimize_angle = false;
blksz = options.render.block_size;

imgpath = [parentpath, filesep, 'stitched'];
maskrpath = [parentpath, filesep, 'masks_region'];
maskwpath = [parentpath, filesep, 'masks_wrinkle'];
ptpairpath = [parentpath, filesep, 'ptpairs'];
outmeshpath = [parentpath, filesep, 'meshes_affine'];
outimgpath = [parentpath, filesep, 'output'];
outxypath = [parentpath, filesep, 'displacement'];
discardpath = [imgpath, filesep, '2nd_pass'];

imglist = dir([imgpath, filesep, '*', ext]);
discardlist = dir([discardpath, filesep, '*', ext]);
discardidx = ismember({imglist.name}, {discardlist.name});
imglist = imglist(~discardidx);


if isempty(imglist)
    return
end
if ~exist(outmeshpath, 'dir')
    mkdir(outmeshpath);
end
if ~exist(outxypath, 'dir')
    mkdir(outxypath);
end
if ~exist(outimgpath, 'dir')
    mkdir(outimgpath);
end
%% filter ptpairs
Nimg = numel(imglist);
matlist = dir([ptpairpath, filesep, '*.mat']);
fixed_sec = false(Nimg, 1);
[Ms_out, cost_history] = optimize_wrapper_multires(imglist, matlist, fixed_sec, options.optimizer, maskwpath, maskrpath);

%%
minxy = nan(Nimg,2);
maxxy = nan(Nimg,2);
for k = 1:Nimg
    minxy(k,:) = min(Ms_out{k}.TR.Points);
    maxxy(k,:) = max(Ms_out{k}.TR.Points);
end
cutratio = 0;
outsz0 = round(fliplr(quantile(maxxy, 1-cutratio, 1) - quantile(minxy, cutratio, 1)));
offst = quantile(minxy, cutratio, 1);
outsz = ceil(outsz0/blksz)*blksz;
offst = offst - fliplr(outsz - outsz0)/2;
% outsz = size(IMG);
% offst = [0,0];
%%
for k = 1:Nimg
    M = Ms_out{k};
    M.TR.Points =  M.TR.Points - offst;
    img = imread([imgpath, filesep, imglist(k).name]);
    if exist([maskrpath, filesep, imglist(k).name], 'file')
        maskr = imread([maskrpath, filesep, imglist(k).name]);
    else
        maskr = ones(size(img), 'uint8');
    end
    if exist([maskwpath, filesep, imglist(k).name], 'file')
        maskw = imread([maskwpath, filesep, imglist(k).name]);
    else
        maskw = zeros(size(img), 'uint8');
    end
    save([outmeshpath, filesep, strrep(imglist(k).name, ext, '.mat')], 'M', 'maskr','maskw', 'outsz');
end

%%
if ~toRender
    return
end
rend_opt = options.render;

%%
pause(5);
block_tform_render_whole_section

% imglist0 = dir('F:\octopus\Sample9-HR-Triad_Tiff\1nm\*.tif');
% outdir0 = 'F:\octopus\Sample9-HR-Triad_Tiff\Sample9-HR-Triad_Aligned\';
% scl = 2;

% for k = 1:Nimg
%     tic
% %     if exist([outimgpath, filesep, imglist(k).name],'file')
% %         continue
% %     end
%     img = imread([imgpath, filesep, imglist(k).name]);
%     load([outmeshpath, filesep, strrep(imglist(k).name, ext, '.mat')])
% 
%     if exist([maskrpath, filesep, imglist(k).name], 'file')
%         maskr = imread([maskrpath, filesep, imglist(k).name]);
%     else
%         maskr = ones(size(img), 'uint8');
%     end
%     if exist([maskwpath, filesep, imglist(k).name], 'file')
%         maskw = uint8(imread([maskwpath, filesep, imglist(k).name]) == 200);
%     else
%         maskw = zeros(size(img), 'uint8');
%     end
%     mxy = elastic_mesh.subregion_movingcoord(M, cat(3,maskr,maskw), rend_opt, [], outsz);
%     % utils.parsave([outxypath, filesep, strrep(imglist(k).name, ext, '.mat')], single(mxy), 'mxy');
%     D = imgaussfilt(elastic_mesh.movcoord2displ(mxy), 15);
%     upscl = 2;
%     imgt = imresize(imwarp(imresize(img, upscl), upscl * imresize(D,upscl), 'FillValues', 255), 1/upscl);
%     imwrite(imgt, [outimgpath, filesep, imglist(k).name]);
% %  
% %     img0 = imread([imglist0(k).folder, filesep, imglist0(k).name]);
% %     D0 = scl * imresize(D, scl);
% %     imgt0 = imwarp(img0, D0);
% %     imwrite(imgt0, [outdir0, filesep, strrep(imglist0(k).name,'.tif','.png')]);
% 
%     toc
% end
