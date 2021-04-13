function mask = recolor_seg(mask, ordered)
    if nargin < 2
        ordered = false;
    end
    if ordered
        [~, ~, idx] = unique(mask(mask>0));
        mask(mask > 0) = idx;
    else
        Nreg = max(mask(:));
        newmap = randperm(255, Nreg);
        mask(mask>0) = newmap(mask(mask>0));
    end
end