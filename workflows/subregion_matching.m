function [ptpairs, cnt, svs, ratio, rgpairs] = subregion_matching(IMG1, IMG2, maskw, maskr, options, A)
    if nargin < 6
        A = [];
    end
    info1 = local_radon.raw_img_to_descriptor_subregion(IMG1, options, {maskw{1}, maskr{1}});
    info2 = local_radon.raw_img_to_descriptor_subregion(IMG2, options, {maskw{2}, maskr{2}});
    matcher_func_name = options.matcher.func;
    matcher_func = str2func(matcher_func_name);
    matcher_params = options.matcher.params;
    ptpairs0 = matcher_func(info1.kps, info2.kps, matcher_params, [], A);
    ptpairs = geometries.filter_matches_subregion(ptpairs0, options.filters);
    [cnt, ratio, svs, rgpairs] = local_paircount(ptpairs0, ptpairs);
end

function [cnt, ratio, svs, rid0u] = local_paircount(ptpairs0, ptpairs)
    rid0 = ptpairs0.region_id;
    rid1 = ptpairs.region_id;
    [rid0u, ~, ic] = unique(rid0, 'rows');
    Np = size(rid0u, 1);
    cnt = nan(Np, 1);
    ratio = nan(Np, 1);
    svs = nan(Np, 1);
    for k = 1:Np
        t0 = sum(ic(:) == k);
        t1 = sum(rid0(:,1) == rid0u(k, 1));
        t2 = sum(rid0(:,2) == rid0u(k, 2));
        ratio(k) = t0 / (min(t1, t2) + 1);
        idxt1 = all(rid1 == rid0u(k,:),2);
        cnt(k) = sum(idxt1);
        if cnt(k) >= 3
            A = geometries.fit_affine(ptpairs.yx1(idxt1,:), ptpairs.yx2(idxt1,:));
            s = svd(A(1:2,1:2));
            svs(k) = max(abs(log(s)));
        end
    end
end