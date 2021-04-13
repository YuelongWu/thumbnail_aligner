function PTPAIRS = stitch_one_section_ptpairs(imglist, YXs, IMGSZs, options, outdir, inverted)
    if nargin < 6 || isempty(inverted)
        inverted = false;
    end
    if nargin < 5
        outdir = [];
    end
    Nstack = 7;
    IMGSTACK = cell(Nstack, 1);
    STACKID = nan(Nstack, 1);
    options = utils.set_default(options, 'min_overlap', 0);
    options = utils.set_default(options, 'scale', 1);

    Nimg = numel(imglist);
    % get a list of overlaps
    if size(IMGSZs, 1) == 1
        IMGSZs = repmat(IMGSZs, Nimg, 1);
    end
    if isempty(IMGSZs)
        IMGSZs = zeros(Nimg, 2);
        for k = 1:Nimg
            if isfield(imglist(k), 'img') && ~isempty(imglist(k).img)
                IMG = imglist(k).img;
            else
                IMG = imread([imglist(k).folder, filesep, imglist(k).name]);
                if inverted
                    IMG = 255 - IMG;
                end
            end
            IMGSZs(k,:) = size(IMG);
        end
    end
    connected = zeros(Nimg, Nimg, 'int8');
    for k = 1:Nimg-1
        bbox0 = local_yx_to_bbox(YXs(k,:), IMGSZs(k,:));
        bbox1 = local_yx_to_bbox(YXs((k+1):end, :), IMGSZs((k+1):end, :));
        [~, overlaps] = rect_intersect(bbox1, bbox0);
        connected(k, (k+1):end) = overlaps > options.min_overlap;
    end
    [kk0, kk1] = find(connected);
    Npp = numel(kk0);
    PTPAIRS = cell(Npp, 1);
    STACKPT = 1;
    for k = 1:Npp
        k0 = kk0(k);
        k1 = kk1(k);

        if ~isempty(outdir)
            outname = [strjoin(erase({imglist(k0).name, imglist(k1).name},'.tif'),'_'),'.mat'];
            if exist([outdir,filesep,outname], 'file')
                continue
            end
        end
        if isfield(imglist(k0), 'img') && ~isempty(imglist(k0).img)
            IMG0 = imglist(k0).img;
        elseif any(STACKID == k0)
            IMG0 = IMGSTACK{find(STACKID == k0,1)};
        else
            IMG0 = imread([imglist(k0).folder, filesep, imglist(k0).name]);
            if inverted
                IMG0 = 255 - IMG0;
            end
            IMGSTACK{STACKPT} = IMG0;
            STACKID(STACKPT) = k0;
            STACKPT = STACKPT + 1;
            if STACKPT > Nstack
                STACKPT = STACKPT - Nstack;
            end
        end
        if isfield(imglist(k1), 'img') && ~isempty(imglist(k1).img)
            IMG1 = imglist(k1).img;
        elseif any(STACKID == k1)
            IMG1 = IMGSTACK{find(STACKID == k1,1)};
        else
            IMG1 = imread([imglist(k1).folder, filesep, imglist(k1).name]);
            if inverted
                IMG1 = 255 - IMG1;
            end
            IMGSTACK{STACKPT} = IMG1;
            STACKID(STACKPT) = k1;
            STACKPT = STACKPT + 1;
            if STACKPT > Nstack
                STACKPT = STACKPT - Nstack;
            end
        end
        YX0 = YXs(k0,:);
        YX1 = YXs(k1,:);
        if options.scale ~= 1
            IMG0 = imresize(IMG0, options.scale);
            IMG1 = imresize(IMG1, options.scale);
            YX0 = YX0 * options.scale;
            YX1 = YX1 * options.scale;
        end
        ptpairs = stitch_two_tiles(IMG0, IMG1, YX0, YX1, options);
        if isempty(ptpairs)
            % warning(['missing overlaps:', imglist(k0).name,'<-->',imglist(k1).name]);
            continue
        end
        if options.scale ~= 1
            ptpairs.yx1 = ptpairs.yx1 / options.scale;
            ptpairs.yx2 = ptpairs.yx2 / options.scale;
        end
        imgnames = {imglist(k0).name, imglist(k1).name};
        strct = struct;
        strct.imgnames = imgnames;
        strct.ptpairs = ptpairs;
        if ~isempty(outdir)
            save([outdir, filesep, outname], 'ptpairs','imgnames');
        end

        PTPAIRS{k} = strct;
    end
    PTPAIRS = vertcat(PTPAIRS{:});
end

function bbox = local_yx_to_bbox(yxs, szs)
    yxs = repelem(yxs, 1, 2);
    Nsz = size(szs, 1);
    szs = [ones(Nsz, 1), szs(:,1), ones(Nsz,1), szs(:,2)];
    bbox = yxs + szs;
end