function ptpairs_flt = filter_matches_subregion(ptpairs, options, L)
    if nargin < 3
        L = [];
    end
    if isfield(options, 'filters')
        options = options.filters;
    end
    if isempty(ptpairs)
        ptpairs_flt = ptpairs;
        return
    end
    options = utils.set_default(options, 'apply_limit', true);
    options = utils.set_default(options, 'rotation_limit', inf);
    options = utils.set_default(options, 'shear_limit', 1);
    options = utils.set_default(options, 'displc_limit', inf);
    filternames = fieldnames(options);
    filteridx = contains(filternames, 'filter');
    filternames = filternames(filteridx);
    [~, ~, ic] = unique(ptpairs(:,'region_id'));
    Ngrp = max(ic(:));
    ptpairs_flt = cell(Ngrp,1);
    for kg = 1:Ngrp
        ptpr = ptpairs(ic == kg, :);
        if size(ptpr,1) < 3
            continue
        end
        for k = 1:numel(filternames)
            fltr = options.(filternames{k});
            fltr_name = fltr.func;
            fltr_func = str2func(fltr_name);
            fltr_params = fltr.params;
            ptpr = fltr_func(ptpr, fltr_params, L);
        end
        Np = size(ptpr, 1);
        if Np <= 2
            continue
        end
        if options.apply_limit
            A = geometries.fit_affine(ptpr.yx1, ptpr.yx2);
            displc = sum((mean(ptpr.yx2, 1) * A(1:2,1:2) + A(3, 1:2) - mean(ptpr.yx2,1)).^2) ^ 0.5;
            if displc > options.displc_limit
                continue
            end
            if ~all(isfinite(A(:)))
                continue;
            end
            [U, S, V] = svd(A(1:2,1:2));
            shear = max(abs(log(diag(S))));
            if shear > options.shear_limit
                continue
            end
            R = U * V';
            theta = acos(R(1));
            if rad2deg(theta) > options.rotation_limit
                continue
            end
        end
        ptpairs_flt{kg} = ptpr;
    end
    ptpairs_flt = cat(1, ptpairs_flt{:});
    if ~isempty(ptpairs_flt)
        ptpairs_flt = geometries.prune_ptpairs(ptpairs_flt);
    end
end