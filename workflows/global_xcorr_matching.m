function ptpairs = global_xcorr_matching(IMG1, IMG2, maskw, maskr, options, A)
    % dyx = yx1 - yx2
    if nargin < 6
        A = eye(3);
    end

    options = utils.set_default(options, 'double_pad', false);
    options = utils.set_default(options, 'apply_mask', false);
    options = utils.set_default(options, 'dis_thresh', 0.5);
    options = utils.set_default(options, 'mask_disk', 1e-3);

    dyx = fliplr(A(3, 1:2));
    if isempty(maskr)
        maskr1 = maskw{1} < 150;
        maskr2 = maskw{2} < 150;
    else
        maskr1 = maskr{1};
        maskr2 = maskr{2};
    end
     if isempty(maskw)
        maskw1 = 255 * (1 - maskr{1});
        maskw2 = 255 * (1 - maskr{2});
    else
        maskw1 = maskw{1};
        maskw2 = maskw{2};
    end
    pre_func_name = options.preprocessing.func;
    pre_func = str2func(pre_func_name);
    pre_params = options.preprocessing.params;
    IMG1f = pre_func(IMG1, pre_params, maskw1);
    IMG2f = pre_func(IMG2, pre_params, maskw2);

    [IMG1f_crp, bbox1] = mask_gen.crop_mask(IMG1f, maskr1);
    [IMG2f_crp, bbox2] = mask_gen.crop_mask(IMG2f, maskr2);

    blksz1 = size(IMG1f_crp);
    blksz2 = size(IMG2f_crp);
    if options.double_pad || any(blksz1 ~= blksz2)
        blksz = max(blksz1, blksz2);
        if options.double_pad
            if max(blksz) / min(blksz) > 3
                [~, paddim] = min(blksz);
                blksz(paddim) = 2 * blksz(paddim);
            else
                paddim = 3;
                blksz = 2 * blksz;
            end
        else
            paddim = 0;
        end
        if any(blksz1 < blksz)
            IMG1f_crp = padarray(IMG1f_crp, blksz - blksz1, 'post');
        end
        if any(blksz2 < blksz)
            IMG2f_crp = padarray(IMG2f_crp, blksz - blksz2, 'post');
        end
    end
    if options.double_pad
        expected_dyxt = [];
    else
        expected_dyxt = -dyx + bbox1(1:2) - bbox2(1:2);
    end
    apply_mask = options.apply_mask;
    mask_disk = options.mask_disk;
    dis_thresh = options.dis_thresh;
    [dyxt , conf] = PMCC.xcorr_fft(IMG1f_crp, IMG2f_crp, mask_disk, expected_dyxt, ...
        [false, false], dis_thresh, apply_mask, paddim);
    dyx_real = -dyxt + bbox1(1:2) -bbox2(1:2);
    bbox = rect_intersect([1,size(IMG1,1), 1, size(IMG1,2)], ...
        repelem(dyx_real,2) + [1, size(IMG2,1), 1,size(IMG2,2)]);
    yx1 = combvec(bbox(1:2), bbox(3:4))';
    myx1 = mean(yx1,1);
    yx1 = 0.5*(yx1 - myx1) + myx1;
    yx2 = yx1 - dyx_real;
    Npt = size(yx1,1);
    conf = conf * ones(Npt,1);
    region_id = ones(Npt, 2);
    rotation = zeros(Npt, 1);
    ptpairs = table(yx1, yx2, region_id, conf, rotation);
 end
