function [ptpairs, Np] = global_matching(IMG1, IMG2, options, varargin)
    info1 = local_radon.raw_img_to_descriptor_global(IMG1, options);
    info2 = local_radon.raw_img_to_descriptor_global(IMG2, options);

    matcher_func_name = options.matcher.func;
    matcher_func = str2func(matcher_func_name);
    matcher_params = options.matcher.params;
    [ptpairs, Mdis, Mrot] = matcher_func(info1.kps, info2.kps, matcher_params);

    match_filters = options.filters;
    filternames = fieldnames(match_filters);
    filteridx = contains(filternames, 'filter');
    filternames = filternames(filteridx);
    for k = 1:numel(filternames)
        fltr = match_filters.(filternames{k});
        fltr_name = fltr.func;
        fltr_func = str2func(fltr_name);
        fltr_params = fltr.params;
        ptpairs = fltr_func(ptpairs, fltr_params);
    end
    Np1 = size(ptpairs, 1);
    xy1 = fliplr(ptpairs.yx1);
    xy2 = fliplr(ptpairs.yx2);
    A = geometries.fit_affine(xy1, xy2);
    ptpairs = matcher_func(info1.kps, info2.kps, matcher_params, [], {Mdis, Mrot, A});
    for k = 1:numel(filternames)
        fltr = match_filters.(filternames{k});
        fltr_name = fltr.func;
        fltr_func = str2func(fltr_name);
        fltr_params = fltr.params;
        ptpairs = fltr_func(ptpairs, fltr_params);
    end
    Np2 = size(ptpairs,1);
    Np = [Np1, Np2];
end