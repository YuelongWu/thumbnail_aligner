addpath('align_config');
parent_path = get_parent_path;
rigidpath = [parent_path, filesep, 'stitched'];
options = parse_yaml(['align_config', filesep, 'wrinkle_mask.yaml']);
ext = '.png';
maskpath = [fileparts(rigidpath), filesep, 'masks_wrinkle'];

if ~exist(maskpath, 'dir')
    mkdir(maskpath)
end
imglist = dir([rigidpath, filesep, '*', ext]);

for k = 1:numel(imglist)
    if exist([maskpath, filesep, imglist(k).name], 'file')
        % continue
    end
    IMG = imread([imglist(k).folder, filesep, imglist(k).name]);
    mask = mask_gen.mask_wrinkle(IMG, options);
    mask(mask <= 200) = 0;
    imwrite(mask, [maskpath, filesep, imglist(k).name]);
end