addpath('align_config');
parent_path = get_parent_path;
wpath = [parent_path, filesep, 'masks_wrinkle'];
options = parse_yaml(['align_config', filesep, 'wrinkle_mask.yaml']);
ext = '.png';
maskpath = [parent_path, filesep, 'masks_region'];

clup = true;
clupsz = 200;

if ~exist(maskpath, 'dir')
    mkdir(maskpath)
end
imglist = dir([wpath, filesep, '*', ext]);

for k = 1:numel(imglist)
    if exist([maskpath, filesep, imglist(k).name], 'file')
        continue
    end
    IMG = imread([imglist(k).folder, filesep, imglist(k).name]);
    if clup
        IMG = mask_gen.fill_small_holes(IMG, clupsz);
    end
    % IMG = imclose(IMG, strel('disk',3));
    mask = mask_gen.partition_section_extend_wrinkle_unique(IMG == 200, IMG ==255, options);
    % mask = mask_gen.partition_section_extend_wrinkle(IMG == 200, IMG ==255, options);
    mask = mask_gen.expand_mask(mask) .* uint8(IMG<255);
    mask = mask_gen.fill_small_holes(mask, clupsz);
    mask = mask .*uint8(IMG<255).*uint8(IMG ~= 200);
    mask = mask_gen.recolor_seg(mask, true);
    imwrite(mask, [maskpath, filesep, imglist(k).name]);
end