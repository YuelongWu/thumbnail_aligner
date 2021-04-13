function [maskw, maskr] = convert_vast_mask(IMGr, IMGw, maskw)
    rids = unique(IMGr(IMGr>0));
    Nrgn = numel(rids);
    maskw = max(single(maskw), 200 * single(IMGw));
    if Nrgn < 2
        maskr = maskw < 255;
        return
    end
    imght = size(IMGr, 1);
    imgwd = size(IMGr, 2);
    Ds = nan(imght, imgwd, Nrgn);
    for k = 1 : Nrgn
        mask0 = IMGr == rids(k);
        Dis1 = bwdistgeodesic(~IMGw, mask0, 'quasi-euclidean');
        Dis2 = bwdist(mask0, 'quasi-euclidean');
        Disd = Dis1 - Dis2;
        Disd = Disd / max(1, max(Disd(:)));
        Ds(:,:,k) = Dis1 .* (Disd + 0.01);
    end
    isolated = all(isnan(Ds), 3);
    [~, maskr] = min(Ds, [], 3);
    maskr(isolated) = 0;
    maskr(maskw > 150) = 0;
    maskr = mask_gen.fill_small_holes(maskr, 200);
end
