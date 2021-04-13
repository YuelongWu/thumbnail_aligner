function maskt = mask_wrinkle(IMG, options)
    %% options
    % wrinkle threshold
    options = utils.set_default(options,'wrk_thresh', [0, 50]);
    options = utils.set_default(options,'dirt_thresh', 10);
    options = utils.set_default(options,'dirt_sz', 5);
    options = utils.set_default(options,'wrk_width', 2);
    options = utils.set_default(options,'wrk_minlen', 100);
    options = utils.set_default(options,'wrk_ratio', 30);
    options = utils.set_default(options,'roi_thresh', [0, inf]);
    options = utils.set_default(options,'roi_erode', 1);
    options = utils.set_default(options,'multicolor', true);

    if numel(options.wrk_thresh) == 1
        options.wrk_thresh = [nan, options.wrk_thresh];
    end
    if numel(options.roi_thresh) == 1
        if options.roi_thresh < 128
            options.roi_thresh = [options.roi_thresh, inf];
        else
            options.roi_thresh = [-inf, options.roi_thresh];
        end
    end

    %%
    mask_roi = (IMG <= options.roi_thresh(1)) | (IMG >= options.roi_thresh(2));
    mask_roi = ~imfill(~mask_roi, 'holes');    
    % mask_roi = mask_roi & (~imclearborder(mask_roi));
    if options.roi_erode > 0
        mask_roi = imdilate(mask_roi, strel('disk', ceil(options.roi_erode)));
    end
    [IMG, bbox] = mask_gen.crop_mask(IMG, ~mask_roi);
    mask = (IMG < min(options.wrk_thresh(2), inf)) & (IMG > max(options.wrk_thresh(1), -inf)) & (~mask_roi(bbox(1):bbox(3), bbox(2):bbox(4)));
    mask = imfill(mask, 'holes');
    if options.wrk_width > 0
        mask = imopen(mask, strel('disk', ceil(options.wrk_width/2)));
    end
    % maskrt = zeros(size(mask));
    % maskln = zeros(size(mask));
    if (options.wrk_minlen > 0) || (options.wrk_ratio > 0)
        CC = bwconncomp(mask, 8);
        Nobj = CC.NumObjects;

        D = bwdist(~mask);
        for k = 1 : Nobj
            idx = CC.PixelIdxList{k};
            wd0 = 2 * mean(D(idx));
            len0 = numel(idx) / wd0;
            ratio0 = len0 / wd0;
            % maskrt(idx) = ratio0;
            % maskln(idx) = len0;
            if (len0 < options.wrk_minlen) || (ratio0 < options.wrk_ratio)
                mask(idx) = false;
            end
        end
    end
    mask_drt = imopen(IMG < options.dirt_thresh, strel('disk', options.dirt_sz));
    mask = uint8(max(100 * mask_drt, 200 * mask));
    maskt = zeros(size(mask_roi), 'uint8');
    maskt(bbox(1):bbox(3), bbox(2):bbox(4)) = mask;
    maskt = max(uint8(255 * mask_roi), maskt);
    if ~options.multicolor
        maskt = maskt == 0;
    end
end