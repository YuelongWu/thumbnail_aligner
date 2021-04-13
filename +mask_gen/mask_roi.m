function mask = mask_roi(IMG, options)
    if isfield(options, 'roi_thresh') && ~isempty(options.roi_thresh)
        threshr = options.roi_thresh;
        if numel(threshr) == 1
            if threshr < 128
                mask = IMG > threshr;
            else
                mask = IMG < threshr;
            end
        else
            mask = (IMG > threshr(1)) & (IMG < threshr(2));
        end
        if isfield(options, 'roi_erode') && (options.roi_erode > 0)
            mask = imerode(mask, strel('disk', round(options.roi_erode)));
        end
    else
        mask = true(size(IMG));
    end
end