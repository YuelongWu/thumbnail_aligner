function M = init_mesh_subregion(M0, mask_region, mask_wrinkle, options)
    if nargin < 4 || isempty(options)
        options = struct;
    end
    options = utils.set_default(options, 'expand_region', true);
    options = utils.set_default(options, 'wrinkle_feather', 50);
    options = utils.set_default(options, 'inregion_soften', true);
    options = utils.set_default(options, 'mass_power', 0);
    options = utils.set_default(options, 'smooth_fun', 'gaussian');
    options = utils.set_default(options, 'stiff_power', 1);
    options = utils.set_default(options, 'cross_region_multiplier', 1);

    if ischar(mask_region)
        mask_region = imread(mask_region);
    end
    mask_region = double(mask_region);
    if ischar(mask_wrinkle)
        mask_wrinkle = imread(mask_wrinkle);
        mask_wrinkle = mask_wrinkle == 200;
    end
    mask_wrinkle = double(mask_wrinkle);
    M = struct;
    M.sideL = M0.sideL;
    M.TR0 = M0.TR;
    yx0 = M0.yx;
    edg0 = M0.edges;
    cList0 = M0.cList;
    M.pt_num0 = size(yx0, 1);
    M.tri_num0 = size(cList0, 1);
    M.edg_num0 = size(edg0, 1);
    region_num = double(max(mask_region(:)));
    if options.mass_power ~= 0
        rmass = histcounts(mask_region(:), 0.5:1:(region_num + 0.5));
        rmass = rmass(:) .^ options.mass_power;
    else
        rmass = ones(region_num, 1);
    end
    if options.params.smooth_gradient
        switch options.smooth_fun
        case 'gaussian'
            M.smooth_kernel = exp(-M0.pt_dis.^2);
        case 'laplacian'
            M.smooth_kernel = exp(-abs(M0.pt_dis));
        otherwise
            M.smooth_kernel = eye(M.pt_num0);
        end
    else
        M.smooth_kernel = [];
    end
    % stacked triangulation
    yx = repmat(yx0, region_num, 1);
    cList = repmat(cList0, region_num, 1);
    indxshift = repmat(0:(region_num-1), M.tri_num0, 1);
    cList = cList + M.pt_num0 * indxshift(:);
    TR = struct;
    TR.Points = fliplr(yx);
    TR.ConnectivityList = cList;
    M.TR = TR;
    M.mass = rmass(:);
    M.pt_num = size(yx, 1);
    M.tri_num = size(cList, 1);
    M.region_num = region_num;
    % vertex id
    yx0_r = min(max(1, round(yx0)), size(mask_region));
    indx = sub2ind(size(mask_region), yx0_r(:,1), yx0_r(:,2));
    if options.expand_region
        mask_region = mask_gen.expand_mask(mask_region);
    end
    M.vtx_id0 = mask_region(indx);

    % edges and stiffness
    if options.wrinkle_feather > 0
        D = bwdist(mask_wrinkle);
        D = D / options.wrinkle_feather;
        mask_wrinkle = min(D, 1).^2;
    else
        mask_wrinkle = 1 - (mask_wrinkle > 0);
    end
    xx = yx0(:, 2);
    yy = yx0(:, 1);
    idx0 = edg0(:, 1);
    idx1 = edg0(:, 2);
    yx0e = [yy(idx0), xx(idx0)];
    yx1e = [yy(idx1), xx(idx1)];
    indx = elastic_mesh.add_lines_indx(size(mask_wrinkle), yx0e, yx1e, true);
    nanidx = all(isnan(indx), 2);
    indx(isnan(indx)) = 1;
    stiffness = min(mask_wrinkle(indx),[],2).^options.stiff_power;
    stiffness(nanidx) = 1;

    [maskg1, maskg2] = elastic_mesh.glue_mask(mask_region);
    g1 = maskg1(indx); g1(nanidx) = 0;
    g2 = maskg2(indx); g2(nanidx) = 0;

    if any(g1(:)) > 0
        g1_flt = g1(g1>0);
        g2_flt = g2(g2>0);
        [m0_flt, ~] = find(g1 > 0);
        tmpmtrx = unique([g1_flt, g2_flt, m0_flt], 'rows');
        g_flt = tmpmtrx(:,1:2);
        glue_idx = tmpmtrx(:,3);
        glue_edg = edg0(glue_idx,:);
        glue_edg = repmat(glue_edg(:),1,2) + repmat(M.pt_num0 * (g_flt - 1),2,1);
        glue_stiffness = options.cross_region_multiplier * repmat(stiffness(glue_idx),2,1);
    else
        glue_edg = [];
        glue_stiffness = [];
    end

    indxshift = repmat(0:(region_num-1), M.edg_num0, 1);
    edg = repmat(edg0, region_num, 1) + M.pt_num0 * indxshift(:);
    edg_id0 = M.vtx_id0(edg0);
    M.edges = [edg; glue_edg];
    if options.inregion_soften
        soft_idx = find(all(edg_id0 > 0, 2) & stiffness < 1 & range(edg_id0, 2) == 0);
        stiff0 = ones(M.edg_num0 * region_num, 1);
        idxt = sub2ind([M.edg_num0, region_num], soft_idx(:), edg_id0(soft_idx(:),1));
        stiff0(idxt) = stiffness(soft_idx);
    else
        stiff0 = ones(M.edg_num0 * region_num,1);
    end
    M.stiffness = [stiff0; glue_stiffness];
    M.L0 = [M.sideL * ones(size(edg,1),1); zeros(size(glue_stiffness))];
end
