function [img, mask0] = pre_imclearborder_log(img0, options, maskt)
    options = utils.set_default(options, 'medfiltsz', 1);
    options = utils.set_default(options, 'feather', 0);
    options = utils.set_default(options, 'feather_small', 0);
    options = utils.set_default(options, 'dilatesz', 0);
    options = utils.set_default(options, 'medfiltsz', 1);
    options = utils.set_default(options, 'sigma', 1);
    options = utils.set_default(options, 'img_power', 1);
    options = utils.set_default(options, 'thresh', 0);
 
    if nargin < 3 || isempty(maskt)
        maskt = 255 * mask_gen.mask_roi(img0, options);
    end
 
    if ischar(maskt)
        maskt = imread(maskt);
    end
    % 0-background 100-artifact, 200-wrinkle, 255-outside
    mask0 = maskt < 150;
    mask0t = ~(maskt == 100);
    if options.medfiltsz > 1
        img0 = medfilt2(img0, [options.medfiltsz, options.medfiltsz], 'symmetric');
    end
    img = imclearborder(img0);
    if options.thresh > 0
        img = max(0, img - options.thresh * max(img(:)));
    end
    img = single(img);
    if options.thresh > 0
        img = max(0, img - options.thresh * options.max_val);
    end
    if options.img_power ~= 1
        img = img.^options.img_power;
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

    if options.sigma > 0
        % img = imgaussfilt(img, options.sigma);
        % img = imfilter(img, [-1, 2, -1], 'replicate') + ...
        %     imfilter(img, [-1; 2; -1], 'replicate');
        img = imfilter(img, fspecial('log', options.sigma*4+1, options.sigma), 'replicate');
    end
end
