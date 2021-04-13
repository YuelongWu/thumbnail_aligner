function [ptpairs, Mdis, Mrot] = descriptor_compare_global(kps1, kps2, options, compare_dim, previous_results)
    options = utils.set_default(options, 'orientation_func', 'local_radon.orient_std_2');
    options = utils.set_default(options, 'compare_method', 'fast');
    options = utils.set_default(options, 'compare_dim', 0); % 0-twoway 1-smaller 2-larger
    options = utils.set_default(options, 'min_corr', 1e-3);
    options = utils.set_default(options, 'refine_max_dis', 100);


    if nargin < 4 || isempty(compare_dim)
        compare_dim = options.compare_dim;
    end
    if nargin > 4 && ~isempty(previous_results)
        Mdis = previous_results{1};
        Mrot = previous_results{2};
        A = previous_results{3};
        D = generate_tformed_distance(kps1, kps2, A);
        Mdis = Mdis .* (D <= options.refine_max_dis);
        previous = true;
    else
        previous = false;
    end
    if ~previous
        if strcmpi(options.compare_method, 'bf')
            [Mdis, Mrot] = distance_brute_force(kps1, kps2);
        elseif strcmpi(options.compare_method, 'fast')
            [Mdis, Mrot] = distance_fast_method(kps1, kps2, options.orientation_func);
        else
            error('Not implemented.')
        end
    end
    Nproj = size(kps1.des,2);
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
    [idx1, idx2] = find(indx);
    yx1 = kps1.yx(idx1, :);
    yx2 = kps2.yx(idx2, :);
    region_id = [kps1.region_id(idx1), kps2.region_id(idx2)];
    conf = ROD(indx) .* Mdis(indx).^2;
    rotation = Mrot(indx) * 360 / Nproj;
    ptpairs = table(yx1, yx2, region_id, conf, rotation);
    ptpairs = sortrows(ptpairs, 'conf', 'descend');
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
