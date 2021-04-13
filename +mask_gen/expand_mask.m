function mask = expand_mask(mask)
    [~, idx] = bwdist(mask > 0);
    mask = mask(idx);
end
