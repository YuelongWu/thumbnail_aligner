addpath('align_config', 'workflows')

old_path = 'F:\octopus\Sample9-HR-Triad_Tiff\verions\version1';
new_path = 'F:\octopus\Sample9-HR-Triad_Tiff\verions\version2';
options = parse_yaml(['align_config', filesep, 'align_same_section.yaml']);
outpath = 'F:\octopus\Sample9-HR-Triad_Tiff\verions\workdir';
ext = '.png';
if ~exist(outpath, 'dir')
    mkdir(outpath)
end

imglist = dir([new_path, filesep, '*', ext]);
if isempty(imglist)
    return
end
% p = gcp;
outimgpath = [outpath, filesep, 'IMG'];
outmeshpath = [outpath, filesep, 'mesh'];
outmxypath = [outpath, filesep, 'mxy0'];
if ~exist(outpath, 'dir')
    mkdir(outpath);
end
if ~exist(outimgpath, 'dir')
    mkdir(outimgpath);
end
if ~exist(outmeshpath, 'dir')
    mkdir(outmeshpath);
end
if ~exist(outmxypath, 'dir')
    mkdir(outmxypath);
end
Nimg = numel(imglist);
for k = 1:Nimg
    try
        if exist([outimgpath, filesep, imglist(k).name], 'file')
            continue
        end
        tic;
        IMG1 = imread([imglist(k).folder, filesep, imglist(k).name]);
        IMG2 = imread([old_path, filesep, imglist(k).name]);
        [IMG2t, M, mxy0] = align_same_section(IMG1, IMG2, options);
        imwrite(imresize(cat(3, IMG2t, IMG1, IMG2t),0.25), [outimgpath, filesep, imglist(k).name]);
        matname = strrep(imglist(k).name, ext, '.mat');
        utils.parsave([outmxypath, filesep, matname], mxy0, 'mxy0');
        utils.parsave([outmeshpath, filesep, matname], M, 'M');
        toc;
    catch ME
        disp(k)
        disp(ME.message)
    end
end
