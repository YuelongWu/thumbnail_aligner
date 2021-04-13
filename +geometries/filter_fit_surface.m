function [ptpairs, dis] = filter_fit_surface(ptpairs, options, L)
    if nargin < 3
        L = [];
    end
    if nargin < 2
        options = [];
    end

    options = utils.set_default(options, 'dis_tol', 5);

    if ~isempty(L) && options.dis_tol < 1
        dis_tol = options.dis_tol * L;
    else
        dis_tol = options.dis_tol;
    end
    yx1 = ptpairs.yx1;
    yx2 = ptpairs.yx2;
    A = geometries.fit_affine(yx1, yx2);
    yx2t = yx2 * A(1:2,1:2) + A(3,1:2);
    dyx = yx2t - yx1;
    if size(ptpairs, 1) < 5
        dis = sum(dyx.^2, 2).^0.5;
    else
        Fx = fit(yx1, dyx(:,2),'lowess');
        Fy = fit(yx1, dyx(:,1),'lowess');
        dxt = Fx(yx1);
        dyt = Fy(yx1);
        dis =((dxt - dyx(:,2)).^2 + (dyt - dyx(:,1)).^2) .^ 0.5;
    end

    idxt = dis <= dis_tol;
    ptpairs = ptpairs(idxt,:);
end
