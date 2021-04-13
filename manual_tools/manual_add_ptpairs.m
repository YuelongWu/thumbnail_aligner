addpath('workflows','align_config')
parentpath = 'D:\forJeff\nervering\TEM_L1_5\Aligned';
algnnames = {'0349.png','0361.png'};

ext = '.png';
imgpath = [parentpath, filesep, 'stitched'];
maskrpath = [parentpath, filesep, 'masks_region'];
maskwpath = [parentpath, filesep, 'masks_wrinkle'];
outpath = [parentpath, filesep, 'ptpairs0'];
imglist = dir([imgpath, filesep, '*', ext]);

listnames = {imglist.name};
kks = nan(size(algnnames));
for k = 1:numel(algnnames)
    kks(k) = find(strcmpi(listnames, algnnames{k}),1);
end

if isempty(imglist)
    return
end

if ~exist(outpath, 'dir')
    mkdir(outpath);
end
warning('off')

for t = 1:size(kks,1)
    kk = kks(t, :);
    k0 = kk(1);
    k1 = kk(2);
    outname = [strrep(imglist(k0).name, ext, ''),'_', strrep(imglist(k1).name, ext, ''),'.mat'];
    if exist([outpath, filesep, outname], 'file')
        continue
    end
    IMG0 = imread([imglist(k0).folder, filesep, imglist(k0).name]);
    maskr0 = imread([maskrpath, filesep, imglist(k0).name]);
    % maskw0 = imread([maskwpath, filesep, imglist(k0).name]);
    IMG1 = imread([imglist(k1).folder, filesep, imglist(k1).name]);
    maskr1 = imread([maskrpath, filesep, imglist(k1).name]);
    % maskw1 = imread([maskwpath, filesep, imglist(k1).name]);
    ptpairs = edit_ptpairs(IMG0, IMG1, [], {maskr0, maskr1});
    imgnames = {imglist(k0).name, imglist(k1).name};
    save([outpath, filesep, outname], 'ptpairs','imgnames');
end