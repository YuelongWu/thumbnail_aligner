function D = movcoord2displ(MC, indxs)
    if nargin > 1 && ~isempty(indxs)
        % indxs: [ymin, ymax, xmin, xmax]
        MC = MC(indxs(1):indxs(2), indxs(3):indxs(4), :);
    end
    imght = size(MC, 1);
    imgwd = size(MC, 2);
    [xx, yy] = meshgrid(1:imgwd, 1:imght);
    MC0 = cat(3,xx,yy);
    D = MC - MC0;
    nanmask = isnan(D);
    if any(nanmask(:))
        [~, idx] = bwdist(~nanmask);
        A = geometries.fit_affine(reshape(D, imght*imgwd, 2), [xx(:), yy(:)], true);
        D0 = [xx(:), yy(:)] * A(1:2, 1:2) + A(3,1:2);
        D0 = reshape(D0, imght, imgwd, 2);
        D1 = D - D0;
        D = D1(idx) + D0; 
    end
end