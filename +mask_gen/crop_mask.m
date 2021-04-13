function [IMGt, bbox] = crop_mask(IMG, mask)
    if all(~mask(:))
        IMGt = [];
        bbox = [];
        return
    end
    m1 = max(mask, [], 1);
    idx1 = find(m1 > 0);
    m2 = max(mask, [], 2);
    idx2 = find(m2 > 0);
    bbox = [idx2(1), idx1(1), idx2(end), idx1(end)];
    IMGt = IMG(bbox(1):bbox(3), bbox(2):bbox(4));
end
