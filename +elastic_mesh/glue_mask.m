function [maskr, maskr_dlt] = glue_mask(maskr, gluesz)
    if nargin < 2
        gluesz = 1;
    end
    maskr_dlt = imdilate(maskr, strel('disk', gluesz));
    mask_brd = maskr ~= maskr_dlt;
    mask_brd = mask_brd & (maskr > 0);
    maskr(~mask_brd) = 0;
    maskr_dlt(~mask_brd) = 0;
end
