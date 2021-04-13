function kps = descriptor_local_radon(img, kps, options)
    persistent rfilters
    if isempty(rfilters)
        rfilters = local_radon.init_radon_filt(options);
    end
    if isfield(rfilters, 'beam_se') && ~isempty(rfilters.beam_se)
        img = imfilter(img, double(rfilters.beam_se), 'replicate');
    end
    if rfilters.filter_blur > 0
        img = imgaussfilt(img, rfilters.filter_blur);
    end
    options = utils.set_default(options, 'regularize_orientaion', false);
    options = utils.set_default(options, 'orientation_func', 'local_radon.orient_2nd_moment');

    yx = kps.yx;
    Npt = size(yx, 1);
    imght = size(img, 1);
    imgwd = size(img, 2);
    proj_num = rfilters.proj_num;
    beam_num = rfilters.beam_num;
    des = zeros(beam_num, Npt, 2 * proj_num, 'single');
    kps_y = yx(:, 1);
    kps_x = yx(:, 2);
    for p = 1:proj_num
        imgp = imfilter(img, rfilters.kernels{p}, 'replicate');
        yb = rfilters.offsets(:, 1, p);
        xb = rfilters.offsets(:, 2, p);
        yy = round(kps_y.' + yb(:));
        xx = round(kps_x.' + xb(:));
        yy = max(1, min(yy, imght));
        xx = max(1, min(xx, imgwd));
        indx = sub2ind(size(imgp), yy, xx);
        des(:,:,p) = imgp(indx);
        yb = rfilters.offsets(:, 1, p + proj_num);
        xb = rfilters.offsets(:, 2, p + proj_num);
        yy = round(kps_y.' + yb(:));
        xx = round(kps_x.' + xb(:));
        yy = max(1, min(yy, imght));
        xx = max(1, min(xx, imgwd));
        indx = sub2ind(size(imgp), yy, xx);
        des(:,:,p + proj_num) = imgp(indx);
    end
    des = permute(des, [1,3,2]); % beam_num x 2proj_num x kp_num
    orientation = [];
    if options.regularize_orientaion
        ort_func = str2func(options.orientation_func);
        [des, orientation] = ort_func(des);
    end
    kps.des = des;
    kps.orientation = orientation(:);
end
