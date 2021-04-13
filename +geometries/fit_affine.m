function [A, R] = fit_affine(XY0, XY1, checknan)
    if nargin < 3
        checknan = false;
    end
    if checknan
        idxnan0 = ~any(isnan(XY0), 2);
        idxnan1 = ~any(isnan(XY1), 2);
        XY0 = XY0(idxnan0 & idxnan1, :);
        XY1 = XY1(idxnan0 & idxnan1, :);
    end
    if isempty(XY0)
        A = eye(3);
        R = A;
        return
    end

    XY0_pad = geometries.padones(XY0);
    XY1_pad = geometries.padones(XY1);
    rnk = min(rank(XY0_pad), rank(XY1_pad));

    if rnk == 1
        A = eye(3);
        A(3,1:2) = XY0 - XY1;
        R = A;
        return
    end
    if rnk == 2
        xyt0 = XY0(1, :) + diff(XY0, 1, 1) * [0, -1; 1,0];
        XY0 = [XY0; xyt0];
        xyt1 = XY1(1, :) + diff(XY1, 1, 1) * [0, -1; 1,0];
        XY1 = [XY1; xyt1];
    end
    warning('off', 'MATLAB:nearlySingularMatrix');
    warning('off', 'MATLAB:singularMatrix');
    warning('off', 'MATLAB:rankDeficientMatrix');
    % center data to increase computational accuracy
    m0 = mean(XY0, 1);
    m1 = mean(XY1, 1);

    XY0 = XY0 - repmat(m0, size(XY0,1), 1);
    XY1 = XY1 - repmat(m1, size(XY1,1), 1);
    XY0_pad = geometries.padones(XY0);
    XY1_pad = geometries.padones(XY1);
    A =  XY1_pad\ XY0_pad;
    if ~all(isfinite(A(:)))
        R = nan(3, 3);
        return
    end
    [U, ~, V] = svd(A(1:2,1:2));
    R = A;
    R(1:2, 1:2) = U * V';
    R(3, 1:2) = R(3, 1:2) + m0 - m1 * R(1:2,1:2);
    R(:, 3) = [0;0;1];
    A(3, 1:2) = A(3, 1:2) + m0 - m1 * A(1:2,1:2);
    A(:, 3) = [0; 0; 1];
end
