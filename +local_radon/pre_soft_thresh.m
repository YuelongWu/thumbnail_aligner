function [img, mask] = pre_soft_thresh(img, options, maskt)
    options = utils.set_default(options, 'max_val', single(intmax(class(img))));
    options = utils.set_default(options, 'thresh', 0);
    options = utils.set_default(options, 'feather', 0);
    options = utils.set_default(options, 'feather_small', 0);
    options = utils.set_default(options, 'dilatesz', 0);
    options = utils.set_default(options, 'medfiltsz', 1);

    if nargin < 3 || isempty(maskt)
        maskt = 255 * (1-mask_gen.mask_roi(img, options));
    end

    if ischar(maskt)
        maskt = imread(maskt);
    end
    % 0-background 100-artifact, 200-wrinkle, 255-outside
    mask0 = maskt < 150;
    mask0t = ~(maskt == 100);

    img = single(img);
    if options.thresh > 0
        img = max(0, img - options.thresh * options.max_val);
    end
    if (options.feather > 0) && (~all(mask0(:)))
        mask0 = bwdist(~mask0) / options.feather;
        mask0 = min(mask0, 1);
    end
    if (options.feather_small > 0) && (~all(mask0t(:)))
        mask0t = bwdist(~mask0t) / options.feather_small;
        mask0t = min(mask0t, 1);
    end
    mask = min(mask0, mask0t);
    img = img .* mask;
    if options.dilatesz > 1
        img = imdilate(img, strel('disk', options.dilatesz));
    end
    if options.medfiltsz > 1
        img = medfilt2(img, [options.medfiltsz, options.medfiltsz], 'symmetric');
    end
end

