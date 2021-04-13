function [IMG2t, M2, mxy2] = align_same_section_PMCC(IMG1, IMG2, options)
    mask1 = true(size(IMG1));
    mask2 = true(size(IMG2));
    fine_options = options.finematch;
    maskr = {ones(size(IMG1)), ones(size(IMG1))};
    maskw = {mask1, mask2};
    % ptpairs0 = [];
    
    C = real(ifft2(fft2(IMG1).*conj(fft2(IMG2))));
    [~, indx] = max(C(:));
    [dy, dx] = ind2sub(size(C), indx);
    dy = dy - 1;
    dx = dx - 1;
    dy = dy - round(dy/size(C,1))*size(C,1);
    dx = dx - round(dx/size(C,2))*size(C,2);
    ptpairs0 = struct;
    ptpairs0.region_id = ones(4,2);
    ptpairs0.yx1 = [0, 0; size(C,1),0; 0, size(C,2); size(C,1), size(C,2)];
    ptpairs0.yx2 = ptpairs0.yx1 - [dy, dx];
    ptpairs = PMCC_matching(IMG1, IMG2, maskw, maskr, fine_options, ptpairs0);
    M2 = [];
    switch options.render.model
        case 'rigid'
            [~, mxy2] = geometries.fit_affine(fliplr(ptpairs.yx1), fliplr(ptpairs.yx2));
            IMG2t = imwarp(IMG2, affine2d(mxy2), 'Outputview', imref2d(size(IMG2)));
            return
        case 'affine'
            [~, mxy2] = geometries.fit_affine(fliplr(ptpairs.yx1), fliplr(ptpairs.yx2));
            IMG2t = imwarp(IMG2, affine2d(mxy2), 'Outputview', imref2d(size(IMG2)));
            imagesc(cat(3, IMG1,IMG2t,IMG1))
            return
        otherwise
            % continue
    end
    R = geometries.fit_affine(fliplr(ptpairs.yx1), fliplr(ptpairs.yx2));
    M0 = elastic_mesh.gen_eqtriang_mesh(size(IMG1), options.optimizer.mesh_space);
    M1 = elastic_mesh.init_mesh_subregion(M0, maskr{1}, maskw{1} == 200, options.optimizer);
    M2 = elastic_mesh.init_mesh_subregion(M0, maskr{2}, maskw{2} == 200, options.optimizer);
    lnk = elastic_mesh.ptpairs_to_links(ptpairs, M1, M2);
    links = elastic_mesh.link_struct(lnk, [1, 2]);
    M2.TR.Points = M2.TR.Points * R(1:2,1:2) + R(3,1:2);
    Ms = {M1, M2};
    Nrep = 20;
    [Ms_out, cost_history0] = elastic_mesh.BGD(Ms, links, options.optimizer.params, [true, false]);
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
    D = imgaussfilt(elastic_mesh.movcoord2displ(mxy2),3);
    IMG2t = imwarp(IMG2, D);
end