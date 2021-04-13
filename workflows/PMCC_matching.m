function [ptpairs, info1, info2] = PMCC_matching(IMG1, IMG2, maskw, maskr, options, ptpairs0)
    info1 = PMCC.raw_img_to_descriptor(IMG1, options, {maskw{1}, maskr{1}});
    info2 = PMCC.raw_img_to_descriptor(IMG2, options, {maskw{2}, maskr{2}});
    matcher_func_name = options.matcher.func;
    matcher_func = str2func(matcher_func_name);
    matcher_params = options.matcher.params;
    L = min(options.detector.params.nb_size);
    if info1.ave_area <= info2.ave_area
        ptpairst = matcher_func(info1.kps, info2.kps, matcher_params, ptpairs0);
        ptpairs = geometries.filter_matches_subregion(ptpairst, options.filters, L);
    else
        ptpairs0 = geometries.flip_ptpairs(ptpairs0);
        ptpairst = matcher_func(info2.kps, info1.kps, matcher_params, ptpairs0);
        ptpairs = geometries.filter_matches_subregion(ptpairst, options.filters, L);
        ptpairs = geometries.flip_ptpairs(ptpairs);
    end
end
