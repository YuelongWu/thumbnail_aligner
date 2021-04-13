function kps = detect_log_blobs(img, options, mask)
    kps = struct;
    options = utils.set_default(options, 'imgpower', 1);
    options = utils.set_default(options, 'medfiltsz', 1);
    options = utils.set_default(options, 'sigma', 3.5);
    options = utils.set_default(options, 'min_dis', ceil(4 * options.sigma));
    options = utils.set_default(options, 'thresh', 0.1);
    options = utils.set_default(options, 'nfeatures', 0);
    options = utils.set_default(options, 'nfeatsubmin', 0.2);
    options = utils.set_default(options, 'nfeatsubgain', 1);

    if (nargin < 3) || isempty(mask)
        mask = 1;
    end
    if ischar(mask)
        mask = imread(mask);
    end
    [rgn_id, ~, ic] = unique(mask(mask > 0));
    N_rgn = numel(rgn_id);
    rgn_areas = histcounts(ic, 0.5:1:(N_rgn+0.5));
    rgn_areas = rgn_areas / sum(rgn_areas(:));
    rgn_areas = max(min(1,  rgn_areas * options.nfeatsubgain), options.nfeatsubmin);

    if options.imgpower ~= 1
        img = img.^(options.imgpower);
    end
    if options.medfiltsz > 1
        img = medfilt2(img, [options.medfiltsz, options.medfiltsz], 'symmetric');
    end

    imgf0 = imgaussfilt(img, options.sigma);
    imgf0 = imfilter(imgf0, [-1, 2, -1], 'replicate') + ...
        imfilter(imgf0, [-1; 2; -1], 'replicate');
    se_r = 50;
    if round(options.min_dis) > 0 && round(options.min_dis) < se_r
        se_r = round(options.min_dis);
    end

    YXs = cell(N_rgn, 1);
    VALs = cell(N_rgn, 1);
    RIDs = cell(N_rgn, 1);
    % idxt = true(N_rgn, 1);
    
    % ridx = 1;
    for kr = 1:N_rgn
        imgf = imgf0 .* (mask == rgn_id(kr));
        if all(imgf(:) <= 0)
            continue
        end
        nfeatures = ceil(options.nfeatures * rgn_areas(kr));
        pkmask = (imgf == imdilate(imgf, strel('disk',se_r))) & ...
            bwareaopen(imgf > (options.thresh * max(imgf(:))), 9);
        [yy, xx] = find(pkmask);
        yx = [yy(:), xx(:)];
        val = imgf(pkmask);
        [val, idx] = sort(val, 'descend', 'MissingPlacement', 'last');
        yx = yx(idx, :);
        Npt = numel(val);
        D = pdist(yx);
        D = D(:);
        toKeep = true(Npt, 1);
        t0 = 0;
        for k = 1:(Npt - 1)
            if toKeep(k)
                tmpD = D((t0 + 1):(t0 + Npt - k));
                toKeep((k+1):end) = toKeep((k+1):end) & (tmpD > options.min_dis);
            end
            t0 = t0 + (Npt - k);
        end
        yx = yx(toKeep, :);
        val = val(toKeep);
        if nfeatures > 0
            if numel(val) > nfeatures
                val = val(1:nfeatures);
                yx = yx(1:nfeatures, :);
            end
        end
        if isempty(yx)
            continue
        end
        YXs{kr} = yx;
        VALs{kr} = val(:);
        RIDs{kr} = single(rgn_id(kr)) * ones(numel(val), 1, 'single');
        % ridx = ridx + 1;
        % if numel(val) == 0
        %     idxt(kr) = false;
        % end
    end
    % YXs = YXs(idxt);
    % VALs = VALs(idxt);
    % RIDs = RIDs(idxt);
    kps.yx = vertcat(YXs{:});
    kps.val = vertcat(VALs{:});
    kps.region_id = vertcat(RIDs{:});
end
