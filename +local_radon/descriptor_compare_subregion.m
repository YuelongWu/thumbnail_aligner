function ptpairs = descriptor_compare_subregion(kps1, kps2, options, compare_dim, A)
    % assume sections has been aligned roughly already
    % kps2 is the moving image (to be broken into pieces)
    options = utils.set_default(options, 'orientation_func', 'local_radon.orient_std_2');
    options = utils.set_default(options, 'compare_method', 'fast');
    options = utils.set_default(options, 'compare_dim', 0); % 0-twoway 1-smaller 2-larger
    options = utils.set_default(options, 'min_corr', 1e-3);
    options = utils.set_default(options, 'expand_whisker', 100);
    options = utils.set_default(options, 'search_whisker', 100);

    if nargin < 4 || isempty(compare_dim)
        compare_dim = options.compare_dim;
    end
    if nargin < 5 || isempty(A)
        A = eye(3);
    end
    idx_pairs = regroup_kps(kps1, kps2, A, options.expand_whisker);
    Nrg = size(idx_pairs, 1);
    YX1 = cell(Nrg, 1);
    YX2 = cell(Nrg, 1);
    RID = cell(Nrg, 1);
    CONF = cell(Nrg, 1);
    ROT = cell(Nrg, 1);
    for kr = 1:Nrg
        kps1t = take_partial_kps(kps1, idx_pairs{kr, 1});
        kps2t = take_partial_kps(kps2, idx_pairs{kr, 2});
        if strcmpi(options.compare_method, 'bf')
            [Mdis, Mrot] = distance_brute_force(kps1t, kps2t);
        elseif strcmpi(options.compare_method, 'fast')
            [Mdis, Mrot] = distance_fast_method(kps1t, kps2t, options.orientation_func);
        else
            error('Not implemented.')
        end
        D = generate_tformed_distance(kps1t, kps2t, A);
        Mdis(D > options.search_whisker) = 0;
        Nproj = size(kps1t.des,2);
        Mdis = max(Mdis, 0.1 * options.min_corr);
        NN = size(Mdis);
        if compare_dim == 0
            ROD1 = compute_ROD(Mdis, 1);
            ROD2 = compute_ROD(Mdis, 2);
            ROD = min(ROD1, ROD2);
        elseif compare_dim == 1
            [~, dim] = min(NN);
            ROD = compute_ROD(Mdis, dim);
        elseif compare_dim == 2
            [~, dim] = max(NN);
            ROD = compute_ROD(Mdis, dim);
        else
            error('Not implemented.')
        end
        indx = (ROD > 1) & (Mdis > options.min_corr);
        if all(~indx)
            continue
        end
        [idx1, idx2] = find(indx);
        YX1{kr} = kps1t.yx(idx1, :);
        YX2{kr} = kps2t.yx(idx2, :);
        RID{kr} = [kps1t.region_id(idx1), kps2t.region_id(idx2)];
        conf = ROD(indx) .* Mdis(indx).^2;
        CONF{kr} = conf(:);
        rotation = Mrot(indx) * 360 / Nproj;
        ROT{kr} = rotation(:);
    end
    yx1 = cat(1, YX1{:});
    yx2 = cat(1, YX2{:});
    region_id = cat(1, RID{:});
    conf = cat(1, CONF{:});
    rotation = cat(1, ROT{:});
    ptpairs = table(yx1, yx2, region_id, conf, rotation);
end

function D = generate_tformed_distance(kps1, kps2, A)
    xy1 = fliplr(kps1.yx);
    xy2 = fliplr(kps2.yx);
    xy2t = xy2 * A(1:2,1:2) + A(3,1:2);
    D = pdist2(xy1, xy2t);
end

function ROD = compute_ROD(Mdis, dim)
    B = min(maxk(Mdis, 2, dim),[],dim);
    ROD = Mdis ./ B;
end

function [Mdis, Mrot] = distance_brute_force(kps1, kps2)
    des1 = kps1.des;
    des2 = kps2.des;
    if isempty(kps1.orientation)
        orient1 = 0;
    else
        orient1 = kps1.orientation(:);
    end
    if isempty(kps2.orientation)
        orient2 = 0;
    else
        orient2 = kps2.orientation(:);
    end
    N1 = size(des1, 3);
    N2 = size(des2, 3);
    Nproj = size(des1, 2);
    Nbeam = size(des1, 1);
    Mdis = -ones(N1, N2, 'single');
    Mrot = zeros(N1, N2, 'single');

    % normalize
    Vdes1 = reshape(des1, Nbeam*Nproj, N1);
    Vdes2 = reshape(des2, Nbeam*Nproj, N2);
    Vdes1 = (Vdes1 - nanmean(Vdes1,1))./nanstd(Vdes1,1,1);
    Vdes2 = (Vdes2 - nanmean(Vdes2,1))./nanstd(Vdes2,1,1);

    for p = 1:Nproj
        if p > 1
            des2 = circshift(des2, 1, 2);
            Vdes2 = reshape(des2, Nbeam*Nproj, N2);
        else
            des2 = reshape(Vdes2, Nbeam, Nproj, N2);
        end
        Mdis1 = Vdes1.' * Vdes2 / (Nproj * Nbeam);
        if p > 1
            Mrot(Mdis1 > Mdis) = (p - 1);
        end
        Mdis = max(Mdis1, Mdis);
    end
    Mrot = Mrot - orient2' + orient1;
    Mrot = Mrot - Nproj*floor(Mrot/Nproj);
end

function [Mdis, Mrot] = distance_fast_method(kps1, kps2, ort_func_name)
    if isempty(kps1.orientation) || isempty(kps2.orientation)
        ort_func = str2func(ort_func_name);
        kps1 = update_orientation(kps1, ort_func);
        kps2 = update_orientation(kps2, ort_func);
    end
    orient1 = kps1.orientation(:);
    orient2 = kps2.orientation(:);
    Mrot = orient1 - orient2';
    % normalize
    N1 = size(kps1.des, 3);
    N2 = size(kps2.des, 3);
    Nproj = size(kps1.des, 2);
    Nbeam = size(kps1.des, 1);
    Vdes1 = reshape(kps1.des, Nbeam*Nproj, N1);
    Vdes2 = reshape(kps2.des, Nbeam*Nproj, N2);
    Vdes1 = (Vdes1 - nanmean(Vdes1,1))./nanstd(Vdes1,1,1);
    Vdes2 = (Vdes2 - nanmean(Vdes2,1))./nanstd(Vdes2,1,1);
    Mdis = Vdes1' * Vdes2 / size(Vdes1, 1);
    Mrot = Mrot - Nproj*floor(Mrot/Nproj);
end

function kps = update_orientation(kps, ort_func)
    if isempty(kps)
        [des, orientation] = ort_func(kps.des);
        kps.des = des;
        kps.orientation = orientation(:);
    end
end

function kps = take_partial_kps(kps, idx)
    if isempty(idx)
        return
    end
    kps.yx = kps.yx(idx, :);
    kps.val = kps.val(idx);
    kps.region_id = kps.region_id(idx);
    if isfield(kps, 'des')
        kps.des = kps.des(:,:,idx);
    end
    if isfield(kps, 'orientation') && ~isempty(kps.orientation)
        kps.orientation = kps.orientation(idx);
    end
end

function idx_pairs = regroup_kps(kps1, kps2, A, expand_whisker)
    region_id2 = kps2.region_id;
    rid2_uniq = unique(region_id2);
    Nid2 = numel(rid2_uniq);
    if Nid2 <= 1
        idx_pairs = {[], []};
        return
    end
    idx_pairs = cell(Nid2, 2);
    D = generate_tformed_distance(kps1, kps2, A);
    idx_mask = ((D == min(D, [], 2)) & (D == min(D, [], 1))) | (D <= expand_whisker);
    for k = 1:Nid2
        idx2 = region_id2 == rid2_uniq(k);
        idx1 = any(idx_mask(:, idx2), 2);
        if all(idx1(:))
            idx1 = [];
        end
        idx_pairs{k, 1} = idx1(:);
        idx_pairs{k, 2} = idx2(:); 
    end
end
