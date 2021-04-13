function [IMG2t, M2, mxy2] = align_same_section(IMG1, IMG2, options)
    % for convert between different versions of alignment
    pre_func1 = str2func(options.prematch.preprocessing.func1);
    param1 = options.prematch.preprocessing.params1;
    [IMG1f, mask1] = pre_func1(IMG1, param1);

    pre_func2 = str2func(options.prematch.preprocessing.func2);
    param2 = options.prematch.preprocessing.params2;
    [IMG2f, mask2] = pre_func2(IMG2, param2);

    detect_func = str2func(options.prematch.detector.func);
    params = options.prematch.detector.params;
    kps1 = detect_func(IMG1f, params);
    kps2 = detect_func(IMG2f, params);

    descriptor_func = str2func(options.prematch.descriptor.func);
    params = options.prematch.descriptor.params;
    kps1 = descriptor_func(IMG1f, kps1, params);
    kps2 = descriptor_func(IMG2f, kps2, params);

    matcher_func = str2func(options.prematch.matcher.func);
    matcher_params = options.prematch.matcher.params;
    ptpairs0 = matcher_func(kps1, kps2, matcher_params);

    match_filters = options.prematch.filters;
    filternames = fieldnames(match_filters);
    filteridx = contains(filternames, 'filter');
    filternames = filternames(filteridx);
    for k = 1:numel(filternames)
        fltr = match_filters.(filternames{k});
        fltr_name = fltr.func;
        fltr_func = str2func(fltr_name);
        fltr_params = fltr.params;
        ptpairs0 = fltr_func(ptpairs0, fltr_params);
    end
    xy1 = fliplr(ptpairs0.yx1);
    xy2 = fliplr(ptpairs0.yx2);
    A = geometries.fit_affine(xy1, xy2);
    % ptpairs0 = matcher_func(kps1, kps2, matcher_params, [], {Mdis, Mrot, A});
    % for k = 1:numel(filternames)
    %     fltr = match_filters.(filternames{k});
    %     fltr_name = fltr.func;
    %     fltr_func = str2func(fltr_name);
    %     fltr_params = fltr.params;
    %     ptpairs0 = fltr_func(ptpairs0, fltr_params);
    % end
    % xy1 = fliplr(ptpairs0.yx1);
    % xy2 = fliplr(ptpairs0.yx2);
    % A = geometries.fit_affine(xy1, xy2);

    IMG2t = imwarp(IMG2, affine2d(A), 'OutputView', imref2d(size(IMG1)));
    xy2t = xy2 * A(1:2, 1:2) + A(3,1:2);
    ptpairs0.yx2 = fliplr(xy2t);
    mask2 = imwarp(mask2, affine2d(A), 'OutputView', imref2d(size(IMG1))) > 0;
    mask2 = 255 - 255 * imdilate(mask2, strel('disk', 2));
    mask1 = 255 - 255 * mask1 > 0.5;

    fine_options = options.finematch;
    maskr = {ones(size(IMG1)), ones(size(IMG1))};
    maskw = {mask1, mask2};
    ptpairs = PMCC_matching(IMG1, IMG2t, maskw, maskr, fine_options, ptpairs0);
    M0 = elastic_mesh.gen_eqtriang_mesh(size(IMG1), options.optimizer.mesh_space);
    M1 = elastic_mesh.init_mesh_subregion(M0, maskr{1}, maskw{1} == 200, options.optimizer);
    M2 = elastic_mesh.init_mesh_subregion(M0, maskr{2}, maskw{2} == 200, options.optimizer);
    lnk = elastic_mesh.ptpairs_to_links(ptpairs, M1, M2);
    links = elastic_mesh.link_struct(lnk, [1, 2]);
    Ms = {M1, M2};
    Nrep = 20;
    [Ms_out, cost_history0] = elastic_mesh.BGD(Ms, links, options.optimizer.params, [true, false], 3);
    cost_history0 = cost_history0(~isnan(cost_history0));
    for k = 1:Nrep
        [Ms_out, cost_history1] = elastic_mesh.BGD(Ms_out, links, options.optimizer.params, [true, false]);
        cost_history1 = cost_history1(~isnan(cost_history1));
        if numel(cost_history1) == numel(cost_history0) && all(cost_history0(:) == cost_history1(:))
            break
        end
        cost_history0 = cost_history1;
    end
    M2 = Ms_out{2};
    mxy2 = elastic_mesh.subregion_movingcoord(M2, maskr{2}, options.render);
    mshp = size(mxy2);
    Ai = inv(A);
    mxy2 = reshape(mxy2, mshp(1) * mshp(2), 2) * Ai(1:2, 1:2) + Ai(3, 1:2);
    mxy2 = reshape(mxy2, mshp);
    D = imgaussfilt(elastic_mesh.movcoord2displ(mxy2),3);
    IMG2t = imwarp(IMG2, D);
end