function [img, mask0] = pre_imreconstruct(img, options, maskt)
    options = utils.set_default(options, 'thresh', 0);
    options = utils.set_default(options, 'feather', 0);
    options = utils.set_default(options, 'mask_dilate', 0);
    options = utils.set_default(options, 'medfiltsz', 1);

    if nargin < 3
        maskt = 255 * (1-mask_gen.mask_roi(img, options));
    end

    if isempty(maskt)
        maskt = 255 * (img == 255);
    end

    if ischar(maskt)
        maskt = imread(maskt);
    end

    % 0-background 100-artifact, 200-wrinkle, 255-outside
    mask0 = uint8(255 * (maskt == 255));

    if options.mask_dilate > 0
        dltsz = ceil(options.mask_dilate);
        mask0 = imdilate(mask0, strel('disk', dltsz));
        mask0(1:dltsz, :) = 255;
        mask0((end-dltsz+1):end, :) = 255;
        mask0(:, 1:dltsz) = 255;
        mask0(:, (end-dltsz+1):end) = 255;
    end

    bck = imreconstruct(mask0, img);
    img = img - bck;
    if options.thresh > 0
        img = max(0, img - options.thresh * max(img(:)));
    end
    img = single(img);
    if (options.feather > 0) && (~all(mask0(:)))
        maskt = bwdist(mask0) / options.feather;
        maskt = min(maskt, 1);
        img = img .* maskt;
    end
    if options.medfiltsz > 1
        img = medfilt2(img, [options.medfiltsz, options.medfiltsz], 'symmetric');
    end
    mask0 = 255 - mask0;
end

