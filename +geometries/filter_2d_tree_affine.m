function ptpairs = filter_2d_tree_affine(ptpairs, options, L)
    if nargin < 3
        L = [];
    end
    options = utils.set_default(options, 'leafsize', 15);
    options = utils.set_default(options, 'dimrange', 250);
    options = utils.set_default(options, 'dis_tol', 20);
    dis_tol = options.dis_tol;
    if dis_tol < 1 && ~isempty(L)
        dis_tol = dis_tol * L;
    end
    IDs = geometries.build_2d_tree(ptpairs.yx1, 1:size(ptpairs,1), options.leafsize, options.dimrange);
    for k = 1:numel(IDs)
        idx = IDs{k};
        idx = idx(:);
        yx1 = ptpairs.yx1(idx, :);
        yx2 = ptpairs.yx2(idx, :);
        A = geometries.fit_affine(yx1, yx2);
        yx2t = yx2 * A(1:2,1:2) + A(3,1:2);
        dis = sum((yx2t - yx1) .^ 2, 2) .^ 0.5;
        idxt = dis < dis_tol;
        if sum(idxt(:)) < 5 && numel(idx) >= 5 
            IDs{k} = [];
        else
            IDs{k} = idx(idxt);
        end
    end
    idx = sort(cat(1, IDs{:}));
    ptpairs = ptpairs(idx, :);
end
