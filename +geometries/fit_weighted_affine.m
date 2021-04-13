function [A, R] = fit_weighted_affine(XY0, XY1, wt, checknan)
    if nargin < 4
        checknan = false;
    end
    if nargin < 3
        wt = ones(size(XY0,1),1);
    end
    wt = wt/sum(wt(:));
    if size(XY0,1) == 0
        A = eye(3);
        R = A;
        return
    end
    if size(XY0,1) == 1
        A = eye(3);
        A(3,1:2) = XY0 - XY1;
        R = A;
        return
    end
    if checknan
        idxnan0 = any(isnan(XY0), 2);
        idxnan1 = any(isnan(XY1), 2);
        XY0 = XY0(idxnan0, :);
        XY1 = XY1(idxnan1, :);
    end

    % center data to increase computational accuracy
    m0 = sum(wt .* XY0, 1);
    m1 = sum(wt .* XY1, 1);

    XY0 = XY0 - repmat(m0, size(XY0,1), 1);
    XY1 = XY1 - repmat(m1, size(XY1,1), 1);
    A = (wt .^ 0.5 .* geometries.padones(XY1)) \ (wt .^ 0.5 .* geometries.padones(XY0));
    A(:, 3) = [0; 0; 1];
    if ~all(isfinite(A(:)))
        R = nan(3, 3);
        return
    end
    if nargout > 1
        [U, ~, V] = svd(A(1:2,1:2));
        R = A;
        R(1:2, 1:2) = U * V';
        R(3, 1:2) = R(3, 1:2) + m0 - m1 * R(1:2,1:2);
        R(:, 3) = [0;0;1];
    end
    A(3, 1:2) = A(3, 1:2) + m0 - m1 * A(1:2,1:2);
    if size(XY0, 1) == 2
        A = R;
    end
end

