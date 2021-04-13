function IMGout = tiles_to_whole_section(run_mode, folderpath, scl)
    if nargin < 3
        scl = 1;
    end
    if nargin < 2
        folderpath = [];
    end
    if nargin < 1
        run_mode = 0;
    end
    IMGout = 'finished.';
    rootpath = '/n/boslfs02/LABS/lichtman_lab/K086';
    switch run_mode
    case 0  % batch folder
        folder_list = dir([rootpath, filesep, 'stitched']);
        outdir = [rootpath, filesep, 'stitched_16nm'];
        idxt = contains({folder_list.name}, '.');
        folder_list = folder_list(~idxt);

        parfor k = 1:numel(folder_list)
            maxNumCompThreads(3);
            try
                outname = [outdir, filesep, folder_list(k).name,'.png'];
                if exist(outname, 'file')
                    continue
                end
                if ~exist([folder_list(k).folder, filesep, folder_list(k).name, filesep, 'rendered'], 'file')
                    disp(['Not finished: ', folder_list(k).name])
                    continue
                end
                IMG = local_tiles_to_whole_section([folder_list(k).folder, filesep, folder_list(k).name], 0.25);
                imwrite(IMG, outname);
            catch ME
                disp(['Error:', folder_list(k).name]);
                disp(ME.message)
            end
        end
    case 1  % downsample
        imglist = dir([rootpath, filesep, 'stitched_16nm', filesep, '*.png']);
        if isempty(imglist)
            return
        end
        outdir1 = [rootpath, filesep, 'stitched_32nm'];
        outdir2 = [rootpath, filesep, 'stitched_128nm'];
        if ~exist(outdir1, 'dir')
            mkdir(outdir1);
        end
        if ~exist(outdir2, 'dir')
            mkdir(outdir2);
        end
        parfor k = 1:numel(imglist)
            outname1 = [outdir1, filesep, imglist(k).name];
            if exist(outname1, 'file')
                continue;
            end
            IMG = imread([imglist(k).folder, filesep, imglist(k).name]);
            IMGt = imresize(IMG, 0.5);
            imwrite(IMGt, [outdir1, filesep, imglist(k).name]);
            IMGtt = imresize(IMGt, 0.25);
            imwrite(IMGtt, [outdir2, filesep, imglist(k).name]);
        end
    otherwise % use as a function
        IMGout = local_tiles_to_whole_section(folderpath, scl);
    end
end

function IMGout = local_tiles_to_whole_section(tilepath, scl, stt)
    if nargin < 3
        stt = [];
    end
    if nargin < 2 || isempty(scl)
        scl = 1;
    end
    ext = '.png';
    imglist = dir([tilepath, filesep, '*', ext]);
    if isempty(imglist)
        disp(['No images found in: ', tilepath]);
        IMGout = [];
        return
    end
    imgnames = erase({imglist.name}, ext);
    outcell = textscan(strjoin(imgnames, newline), '%*s%f%f', 'Delimiter', {'_tr','_tc'});
    r0 = outcell{1};
    c0 = outcell{2};

    if isempty(stt)
        stt_0 = min(min(r0(:)), min(c0(:)));
        r0 = r0 - stt_0 + 1;
        c0 = c0 - stt_0 + 1;
        Nr = max(r0(:));
        Nc = max(c0(:));
    elseif numel(stt) == 4
        % [rmin, rmax, cmin, cmax]
        r0 = r0 - stt(1) + 1;
        c0 = c0 - stt(3) + 1;
        idxt = (r0 > 0) & (c0 > 0) & (r0 <= stt(2) - stt(1) + 1) & (c0 <= stt(4) - stt(3) + 1);
        r0 = r0(idxt);
        c0 = c0(idxt);
        imglist = imglist(idxt);
        Nr = max(r0(:));
        Nc = max(c0(:));
    else
        stt_0 = stt(1);
        r0 = r0 - stt_0 + 1;
        c0 = c0 - stt_0 + 1;
        Nr = max(r0(:));
        Nc = max(c0(:));
        idxt = (r0 > 0) & (c0 > 0);
        r0 = r0(idxt);
        c0 = c0(idxt);
        imglist = imglist(idxt);
    end
    [~, idxt] = sort(r0 + c0);
    r0 = r0(idxt);
    c0 = c0(idxt);
    imglist = imglist(idxt);

    IMG = imread([imglist(1).folder, filesep, imglist(1).name]);
    if scl ~= 1
        IMG = imresize(IMG, scl);
    end
    imght = size(IMG, 1);
    imgwd = size(IMG, 2);
    outsz = [size(IMG, 1) * Nr, size(IMG,2) * Nc, size(IMG, 3)];
    IMGout = 255 * ones(outsz, class(IMG));

    yy0 = (r0(1) - 1) * imght + 1;
    yy1 = yy0 + size(IMG, 1) - 1;
    xx0 = (c0(1) - 1) * imgwd + 1;
    xx1 = xx0 + size(IMG, 2) - 1;
    IMGout(yy0:yy1, xx0:xx1,:) = IMG;

    for k = 2:numel(imglist)
        IMG = imread([imglist(k).folder, filesep, imglist(k).name]);
        if scl ~= 1
            IMG = imresize(IMG, scl);
        end
        yy0 = (r0(k) - 1) * imght + 1;
        yy1 = yy0 + size(IMG, 1) - 1;
        xx0 = (c0(k) - 1) * imgwd + 1;
        xx1 = xx0 + size(IMG, 2) - 1;
        IMGout(yy0:yy1, xx0:xx1,:) = IMG;
    end
end
