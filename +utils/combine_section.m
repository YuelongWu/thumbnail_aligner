function IMG = combine_section(imgpath, expr, bbox, rowfirst, tilesz, tilenum, sttpt)
    % bbox: [ymin, ymax, xmin, xmax]
    if nargin < 7
        sttpt = 1;
    end
    if nargin < 6
        tilenum = [];
    end
    if nargin < 5
        tilesz = [];
    end
    if nargin < 4
        rowfirst = true;
    end
    if nargin < 3
        bbox = [];
    end
    imglist = dir(imgpath);
    if isempty(imglist)
        IMG = 0;
        return
    end
    IMG0 = imread([imglist(1).folder, filesep, imglist(1).name]);
    tp = class(IMG0);
    imgnames = {imglist.name};
    rc = utils.parse_row_col(imgnames(:), expr, rowfirst, sttpt);
    if isempty(tilenum)
        Nr0 = max(rc(:,1));
        Nc0 = max(rc(:,2));
    else
        Nr0 = tilenum(1);
        Nc0 = tilenum(2);
    end
   
    if isempty(tilesz)
        [~, idxt] = min(sum(rc, 2)); % tile with minimum indices
        tile = imread([imglist(idxt).folder, filesep, imglist(idxt).name]);
        imght = size(tile, 1);
        imgwd = size(tile, 2);
    else
        imght = tilesz(1);
        imgwd = tilesz(2);
    end

    if ~isempty(bbox)
        minr = max(1, floor((bbox(1) - 1) / imght) + 1);
        minc = max(1, floor((bbox(3) - 1) / imgwd) + 1);
        maxr = min(Nr0, ceil(bbox(2) / imght));
        maxc = min(Nc0, ceil(bbox(4) / imgwd));
    else
        minr = 1;
        minc = 1;
        maxr = Nr0;
        maxc = Nc0;
    end
    Nr = maxr - minr + 1;
    Nc = maxc - minc + 1;
    tilecell = cell(Nr, Nc);
    missingnum = 0;
    for r = 1:Nr
        for c = 1:Nc
            r0 = r + minr - 1;
            c0 = c + minc - 1;
            idxt = find(all(rc == [r0, c0], 2), 1);
            if isempty(idxt)
                tilecell{r,c} = zeros(imght, imgwd, tp);
                missingnum = missingnum + 1;
            else
                tilecell{r,c} = imread([imglist(idxt).folder, filesep, imglist(idxt).name]);
            end
        end
    end
    if missingnum > 0
        fprintf('missing %d tiles | ', missingnum);
    end
    IMG = cell2mat(tilecell);
    if ~isempty(bbox)
        bbox([1, 2]) = bbox([1, 2]) - (minr - 1) * imght;
        bbox([3, 4]) = bbox([3, 4]) - (minc - 1) * imgwd;
        bbox([2, 4]) = min(bbox([2, 4]), [size(IMG,1), size(IMG,2)]);
        IMG = IMG(bbox(1):bbox(2), bbox(3):bbox(4));
    end
end
    