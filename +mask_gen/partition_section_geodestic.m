function CC = partition_section_geodestic(mask, seed, diff_thresh, area_thresh)
    scl = 0.5;
    inline_thresh = 3;
    if nargin < 2
        seed = init_seed_mask(mask, scl);
    end
    if numel(seed) == 2
        idx = sub2ind(size(mask), seed(1), seed(2));
        seed = false(size(mask));
        seed(idx) = true;
    end
    if nargin < 3
        diff_thresh = 50;
    end
    if nargin < 4
        area_thresh = 1000;
    end
    
    WT = cell(100, 1);

    maskt = imresize(imresize(single(mask), scl, 'bilinear') > 0, size(mask), 'nearest');
    maskt = imerode(maskt, ones(3));
    [~, idxt] = bwdist(maskt);
    covered_dis = nan(size(mask));
    covered_dis(maskt) = -inf;
    t0 = 1;
    while nansum(covered_dis(:) < -diff_thresh) > area_thresh
        D1 = mask_gen.geodesic_dis_diff(mask, seed, scl);
        D1(~maskt) = nan;
        T = D1 < inline_thresh;
        T1 = bwdist(~T);
        D1 = mask_gen.geodesic_dis_diff(mask, T1 >= 0.5*max(T1(:)), scl);
        T = D1 < inline_thresh;
        % T = T(idxt);
        D2 = bwdist(~T, 'quasi-euclidean');
        D = D2 - D1;
        covered_dis = max(covered_dis, D);
        [maxval, idxm] = min(covered_dis(:));
        seed = false(size(mask));
        if isinf(maxval)
            masktt = imfill(~maskt, idxm) & maskt;
            Dt = bwdist(~masktt);
            [~, idxm] = max(Dt(:));
        end
        seed(idxm) = true;
        WT{t0} = D;
        t0 = t0 + 1;
    end
    WT = cat(3, WT{:});
    [~, CC] = max(WT, [], 3);
    CC = CC(idxt) .* mask;
    CC = mask_gen.fill_small_holes(CC, area_thresh);
end

function seed = init_seed_mask(mask, scl)
    vidx = find(mask(:));
    Npix = numel(vidx);
    idx = randi(Npix);
    rand_seed = false(size(mask));
    rand_seed(vidx(idx)) = true;
    D = mask_gen.geodesic_dis_diff(mask, rand_seed, scl);
    D(~mask) = nan;
    seed = false(size(mask));
    [maxval, idx] = max(D(:));

    if isinf(maxval)
        maskt = imfill(~mask, idx) & mask;
        D = bwdist(~maskt);
        [~, idx] = max(D(:));
    end
    seed(idx) = true;
end
