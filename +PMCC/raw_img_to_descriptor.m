function out = raw_img_to_descriptor(IMG, options, masks)
    if ~isstruct(IMG)
        out = struct;
        out.img = IMG;
    else
        out = IMG;
    end
    if ~isfield(out, 'imgf')
        mask0 = masks{1};
        out.maskw = mask0;
        pre_func_name = options.preprocessing.func;
        pre_func = str2func(pre_func_name);
        pre_params = options.preprocessing.params;
        out.imgf = pre_func(out.img, pre_params, mask0);
    end
    if ~isfield(out, 'kps')
        mask1 = masks{2};
        if ischar(mask1)
            mask1 = imread(mask1);
        end
        mask1 = single(mask1);
        rarea = histcounts(mask1(:), 0.5:1:(max(mask1(:))+0.5));
        out.maskr = mask1;
        out.ave_area = mean(rarea(rarea > 0));
        detect_func_name = options.detector.func;
        detect_func = str2func(detect_func_name);
        detect_params = options.detector.params;
        out.kps = detect_func(out.imgf, detect_params, mask1);
    end
end