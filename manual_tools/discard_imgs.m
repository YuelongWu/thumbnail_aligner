function discard_imgs(imgpath, restart, tgtname)
    addpath('align_config');
    if nargin < 3
        tgtname = '2nd_round';
    end
    if nargin < 2
        restart = false;
    end
    if nargin < 1 || isempty(imgpath)
        parent_path = get_parent_path;
        imgpath = [parent_path, filesep, 'stitched'];
    end
    discardpath = [imgpath, filesep, tgtname];
    imglist = dir([imgpath, filesep, '*.png']);

    if ~exist(discardpath, 'dir')
        mkdir(discardpath);
    end
    Nimg = length(imglist);
    crnt_idx = 1;
    dcrdlist = dir([discardpath, filesep, '*.png']);
    if ~restart && ~isempty(dcrdlist)
        Dname = {dcrdlist.name};
        IMGname = {imglist.name};
        Tname = [Dname(:); IMGname(:)];
        Tname = unique(Tname);
        Tname = sort(Tname);
        [~, ib] = ismember(Dname,Tname);
        ib = max(ib(:));
        if ib == numel(Tname)
            disp('previously done.')
            return
        end
        if ib > 0
            tname = Tname{ib+1};
            idxt = strcmpi(IMGname, tname);
            crnt_idx = find(idxt,1);
            if isempty(crnt_idx)
                crnt_idx = 1;
            end
        end
    end
    hfig = figure('Name', 'Select key frames (z:previous x:next)', ...
        'NumberTitle', 'off', 'Unit','normalized','Position',[0.1,0.1,0.8,0.8]);
    colormap(gray);
    ax = axes(hfig);
    while crnt_idx <= Nimg
        localShowIMG(imglist, crnt_idx, ax);
        w = waitforbuttonpress;
        if w
            switch lower(get(hfig,'CurrentCharacter'))
            case {'x', char(29), char(13)}
                crnt_idx = crnt_idx + 1;
                if crnt_idx > Nimg
                    return
                end
            case {'z', char(28)}
                crnt_idx = crnt_idx - 1;
                if crnt_idx < 1
                    crnt_idx = 1;
                end
            case 'q'
                close(hfig)
                return
            case 'd'
                movefile([imglist(crnt_idx).folder, filesep, imglist(crnt_idx).name],...
                    [discardpath, filesep, imglist(crnt_idx).name]);
                imglist = dir([imgpath, filesep, '*.png']);
                Nimg = numel(imglist);
           otherwise
                continue
            end
        end
    end
end

function IMG = localShowIMG(imglist, crnt_idx, ax)
    if nargin < 4
        ax = gca;
    end
    k = crnt_idx;
    IMG = imread([imglist(k).folder, filesep, imglist(k).name]);
    IMGt = imresize(IMG, 1);
    imagesc(ax, IMGt);
    axis(ax, 'off');
    axis(ax, 'equal');
    title([strrep(imglist(k).name, '_', '\_'), ...
        ' (progress: ',num2str(crnt_idx),'/',num2str(numel(imglist)),')']);
end