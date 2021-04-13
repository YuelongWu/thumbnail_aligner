function ptpairs = filter_global_ransac(ptpairs, options, L)
    if nargin < 3
        L = [];
    end
    options = utils.set_default(options, 'max_iter', inf);
    options = utils.set_default(options, 'dis_tol', 10);
    options = utils.set_default(options, 'match_exit_num', inf);
    options = utils.set_default(options, 'match_exit_ratio', 1);
    dis_tol = options.dis_tol;
    if dis_tol < 1 && ~isempty(L)
        dis_tol = dis_tol * L;
    end
    exhaust_ratio = 0.5; % part of the time spent on exhaustive searching
    Np = size(ptpairs, 1);
    if Np < 4
        return
    end
    Ncmb = nchoosek(Np, 3);
    exhaust_ratio = max(exhaust_ratio, min(1, options.max_iter/Ncmb));
    max_iter = min(Ncmb, options.max_iter);
    match_exit_num = min(options.match_exit_num, Np);
    exhaust_quota = max_iter * exhaust_ratio;
    N0 = 4;
    for k = 4:Np
        if nchoosek(k, 3) <= exhaust_quota
            N0 = k;
        else
            break
        end
    end
    Nexh = nchoosek(N0, 3);
    combs = nchoosek(1:N0, 3);
    yx1 = ptpairs.yx1;
    yx2 = ptpairs.yx2;
    hit = false;
    countdown = 20;
    inlier_num = 0;
    inlier_idx = false(Np, 1);
    for  k = 1:max_iter
        if k <= Nexh
            comb = combs(k, :);
        else
            comb = randperm(Np, 3);
            if all(comb <= N0)
                comb(1) = N0 + randi(Np - N0);
            end
        end
        yx1sp = yx1(comb, :);
        yx2sp = yx2(comb, :);
        if numel(unique(yx1sp, 'rows')) ~= numel(yx1sp) ||  numel(unique(yx2sp, 'rows')) ~= numel(yx2sp)
            continue
        end
        A = geometries.fit_affine(yx1sp, yx2sp);
        if ~all(isfinite(A(:)))
            continue
        end
        dyx = yx2 * A(1:2, 1:2) + A(3, 1:2) - yx1;
        inliers = (sum(dyx.^2, 2)).^0.5 <= dis_tol;
        inlier_cnt = sum(inliers);
        if inlier_cnt < 3
            continue
        end
        svs = svd(A(1:2,1:2));
        svms = min(exp(-abs(log(svs))));
        if inlier_cnt * svms > inlier_num
            inlier_num = (inlier_cnt) * svms;
            inlier_idx = inliers;
            if ~hit && inlier_cnt >= max(match_exit_num, Np * options.match_exit_ratio)
                hit = true;
            end
        end
        if hit
            countdown = countdown - 1;
        end
        if countdown == 0
            break
        end
    end
    A = geometries.fit_affine(yx1(inlier_idx,:), yx2(inlier_idx,:));
    dyx = yx2 * A(1:2, 1:2) + A(3, 1:2) - yx1;
    inliers = (sum(dyx.^2, 2)).^0.5 <= dis_tol;
    ptpairs = ptpairs(inliers, :);
end
