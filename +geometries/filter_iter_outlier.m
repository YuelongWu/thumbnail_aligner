function ptpairs = filter_iter_outlier(ptpairs, options, varargin)
    if nargin < 2
        options = [];
    end
    options = utils.set_default(options, 'outlier_thresh', 5);
    options = utils.set_default(options, 'max_delete', 5);
    options = utils.set_default(options, 'delete_batch', 1);
    options = utils.set_default(options, 'min_remain', 3);
    options = utils.set_default(options, 'std_relax', 0.5);
    options = utils.set_default(options, 'upper_shoulder', 15);
    yx1 = ptpairs.yx1;
    yx2 = ptpairs.yx2;
    idxt = true(size(yx1,1),1);
    while (sum(~idxt) < options.max_delete) && sum(idxt) > options.min_remain
        A = geometries.fit_affine(yx1(idxt,:), yx2(idxt,:));
        yx2t = yx2 * A(1:2,1:2) + A(3,1:2);
        dis = sum((yx2t - yx1).^2,2).^0.5;
        if isempty(options.upper_shoulder) || isinf(options.upper_shoulder)
            m = mean(dis(idxt));
        else
            m = min(maxk(dis(idxt), options.upper_shoulder));
        end
        cc = (dis - m).*(idxt) / (std(dis(idxt)) + options.std_relax);
        if all(cc < options.outlier_thresh)
            break
        end
        throwaway = (cc >= options.outlier_thresh) & ...
            (cc >= min(maxk(cc,options.delete_batch)));
        idxt(throwaway) = false;
    end
    ptpairs = ptpairs(idxt, :);
end
