function [region_id, A, R] = prematch2affine(ptpairs0)
    rids = ptpairs0.region_id;
    [region_id, ~, ic] = unique(rids, 'rows');
    Nrp = size(region_id, 1);
    idxt = true(Nrp, 1);
    A = repmat(eye(3), 1, 1, Nrp);
    R = A;
    for k = 1:Nrp
        yx1 = ptpairs0.yx1(ic == k,:);
        yx2 = ptpairs0.yx2(ic == k,:);
        [A0, R0] = geometries.fit_affine(yx2, yx1);
        if ~all(isfinite(A0(:)))
            idxt(k) = false;
        else
            A(:,:,k) = A0;
            R(:,:,k) = R0;
        end
    end
    region_id = region_id(idxt, :);
    A = A(:,:,idxt);
    R = R(:,:,idxt);
end
