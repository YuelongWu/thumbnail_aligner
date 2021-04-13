function [img, mask0] = pre_imclearborder(img0, options)
    options = utils.set_default(options, 'medfiltsz', 1);
    mask0 = img0 > 250;
    mask0 = mask0 & (~imclearborder(mask0));
    mask0 = ~mask0;
    img = imclearborder(img0);
    if options.thresh > 0
        img = max(0, img - options.thresh * max(img(:)));
    end
    img = single(img);
    if options.medfiltsz > 1
        img = medfilt2(img, [options.medfiltsz, options.medfiltsz], 'symmetric');
    end
end


