function out = raw_img_to_descriptor_subregion(IMG, options, mask_dirs)
    if ~isstruct(IMG)
        out = struct;
        out.img = IMG;
    else
        out = IMG;
    end
    if ~isfield(out, 'imgf')
        mask0 = mask_dirs{1};
        out.maskw = mask0;
        pre_func_name = options.preprocessing.func;
        pre_func = str2func(pre_func_name);
        pre_params = options.preprocessing.params;
        out.imgf = pre_func(out.img, pre_params, mask0);
    end
    if ~isfield(out, 'kps')
        mask1 = mask_dirs{2};
        out.maskr = mask1;
        detect_func_name = options.detector.func;
        detect_func = str2func(detect_func_name);
        detect_params = options.detector.params;
        out.kps = detect_func(out.imgf, detect_params, mask1);
    end
    if ~isfield(out.kps, 'des')
        descriptor_func_name = options.descriptor.func;
        descriptor_func = str2func(descriptor_func_name);
        descriptor_params = options.descriptor.params;
        out.kps = descriptor_func(out.imgf, out.kps, descriptor_params);
    end
end