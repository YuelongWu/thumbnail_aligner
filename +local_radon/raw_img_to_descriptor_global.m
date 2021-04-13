function out = raw_img_to_descriptor_global(IMG, options)
    if ~isstruct(IMG)
        out = struct;
        out.img = IMG;
    else
        out = IMG;
    end
    if ~isfield(out, 'imgf')
        mask_func_name = options.mask_gen.func;
        mask_func = str2func(mask_func_name);
        mask_params = options.mask_gen.params;
        mask = mask_func(out.img, mask_params);

        pre_func_name = options.preprocessing.func;
        pre_func = str2func(pre_func_name);
        pre_params = options.preprocessing.params;
        out.imgf = pre_func(out.img, pre_params, mask);
    end
    if ~isfield(out, 'kps')
        detect_func_name = options.detector.func;
        detect_func = str2func(detect_func_name);
        detect_params = options.detector.params;
        out.kps = detect_func(out.imgf, detect_params);
    end
    if ~isfield(out.kps, 'des')
        descriptor_func_name = options.descriptor.func;
        descriptor_func = str2func(descriptor_func_name);
        descriptor_params = options.descriptor.params;
        out.kps = descriptor_func(out.imgf, out.kps, descriptor_params);
    end
end