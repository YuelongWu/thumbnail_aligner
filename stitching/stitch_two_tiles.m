function ptpairs = stitch_two_tiles(IMG1, IMG2, YX1, YX2, options)
    [mask1, mask2] = gen_overlap_mask(YX1, YX2, size(IMG1), size(IMG2), options.mask_gen);
    A = eye(3);
    A(3, 1:2) = fliplr(YX2 - YX1);
    mask1w = 255*(1-mask1);
    mask2w = 255*(1-mask2);

    % hydra specific
    m1 = imdilate(IMG1 == 0, ones(5));
    m2 = imdilate(IMG2 == 0, ones(5));
    mask1w(m1) = max(100, mask1w(m1));
    mask2w(m2) = max(100, mask2w(m2));

    maskw = {mask1w, mask2w};
    maskr = {mask1, mask2};
    yx2 = combvec([1, size(mask2,1)], [1, size(mask2,2)])';
    ptpairs = local_create_dummy_ptpair(yx2, A);

    if isfield(options, 'prematch')
        prem_func_name = options.prematch.func;
        prem_func = str2func(prem_func_name);
        prem_params = options.prematch.params;
        ptpairs = prem_func(IMG1, IMG2, maskw, maskr, prem_params, A);
        ptpairs = ptpairs(ptpairs.conf > options.prematch.min_conf, :);
        if isempty(ptpairs)
            return
        end
    end

    if isfield(options, 'finematch')
        fine_func_name = options.finematch.func;
        fine_func = str2func(fine_func_name);
        fine_params = options.finematch.params;
        ptpairs = fine_func(IMG1, IMG2, maskw, maskr, fine_params, ptpairs);
    end
end

function ptpairs = local_create_dummy_ptpair(yx2, A)
    yx1 = yx2 * A(2:-1:1,2:-1:1) + A(3, 2:-1:1);
    Npt = size(yx1,1);
    conf =  zeros(Npt,1);
    region_id = ones(Npt, 2);
    rotation = zeros(Npt, 1);
    ptpairs = table(yx1, yx2, region_id, conf, rotation);
end
