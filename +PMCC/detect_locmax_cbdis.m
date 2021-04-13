function kps = detect_locmax_cbdis(img0, options, mask)
    kps = struct;
    options = utils.set_default(options, 'nfeature', 5000);
    options = utils.set_default(options, 'min_dis', 100);
    options = utils.set_default(options, 'thresh', 0.05);
    options = utils.set_default(options, 'nfeatsubmin', 0.2);
    options = utils.set_default(options, 'nfeatsubgain', 1.25);
    options = utils.set_default(options, 'nb_size', 201);
    options = utils.set_default(options, 'fft', false);
    options = utils.set_default(options, 'use_abs', true);
    options = utils.set_default(options, 'pad_des', false);

    if (nargin < 3) || isempty(mask)
        mask = 1;
    end
    if ischar(mask)
        mask = imread(mask);
    end
    mask = cast(mask, class(img0));
    [rgn_id, ~, ic] = unique(mask(mask > 0));
    N_rgn = numel(rgn_id);
    rgn_areas = histcounts(ic, 0.5:1:(N_rgn+0.5));
    rgn_ratio = rgn_areas / sum(rgn_areas(:));
    rgn_ratio = max(min(1,  rgn_ratio * options.nfeatsubgain), options.nfeatsubmin);

    min_dis = options.min_dis;
    nb_size = options.nb_size;
    [~, ia_g] = unique([min_dis(:), nb_size(:)], 'rows');
    min_dis = min_dis(ia_g);
    nb_size = nb_size(ia_g);

    kps.min_dis = min_dis;
    kps.nb_size = nb_size;

    nfeatures = options.nfeature;
    if numel(nfeatures) < numel(min_dis)
        nfeatures = nfeatures(1) * ones(size(min_dis));
    end
    N_dis = numel(min_dis);
    YXs = cell(N_rgn, N_dis);
    VALs = cell(N_rgn, N_dis);
    RIDs = cell(N_rgn, N_dis);
    DESs = cell(N_rgn, N_dis);

    imght = size(img0,1);
    imgwd = size(img0,2);
    if options.use_abs
        imgp0 = abs(img0);
    else
        imgp0 = img0;
    end
    for kr = 1:N_rgn
        img = img0 .* (mask == rgn_id(kr));
        imgp = imgp0 .* (mask == rgn_id(kr));
        if all(img(:) <= 0)
            continue
        end
        nfeat = ceil(nfeatures * rgn_ratio(kr));
        for kd = 1:N_dis
            imgf = imdilate(imdilate(imgp, ones(min_dis(kd),1)), ones(1, min_dis(kd)));
            pkmask = imgp == imgf & imgp > (options.thresh * max(imgp(:)));
            [yy, xx] = find(pkmask);
            yx = [yy(:), xx(:)];
            val = imgf(pkmask);
            [val, idx] = sort(val, 'descend', 'MissingPlacement', 'last','ComparisonMethod','abs');
            yx = yx(idx, :);
            if nfeat > 0
                if numel(val) > nfeat
                    val = val(1:nfeat);
                    yx = yx(1:nfeat, :);
                end
            end
            if isempty(yx)
                continue
            end
            YXs{kr, kd} = yx;
            VALs{kr, kd} = val(:);
            RIDs{kr, kd} = rgn_id(kr) * ones(numel(val), 1);

            yy = permute(yx(:,1), [2,3,1]);
            xx = permute(yx(:,2), [2,3,1]);
            wsz = round(nb_size(kd) / 2);
            [xg, yg] = meshgrid(-wsz:1:wsz, -wsz:1:wsz);
            xx = xx + xg;
            yy = yy + yg;
            xx = min(max(xx, 1),imgwd);
            yy = min(max(yy, 1),imght);
            indx = sub2ind([imght, imgwd], yy, xx);
            des0 = img(indx);
            if options.pad_des
                des0 = padarray(des0, floor([size(des0,1), size(des0,2), 0]/2), 'both');
            end
            DESs{kr, kd} = des0;
        end
    end
    kps.yx = cell(N_dis,1);
    kps.val = cell(N_dis,1);
    kps.region_id = cell(N_dis,1);
    kps.des = cell(N_dis,1);
    if numel(options.fft) == 1
        kps.ffted = repmat(options.fft, N_dis, 1);
    else
        kps.ffted = options.fft(ia_g);
    end
    for kd = 1:N_dis
        kps.yx{kd} = vertcat(YXs{:, kd});
        kps.val{kd} = vertcat(VALs{:, kd});
        kps.region_id{kd} = vertcat(RIDs{:,kd});
        if kps.ffted(kd)
            kps.des{kd} = fft2(cat(3, DESs{:, kd}));
        else
            kps.des{kd} = cat(3, DESs{:, kd});
        end
    end
end
