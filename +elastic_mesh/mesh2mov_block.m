function MC = mesh2mov_block(M, masks, options, rid, outsz, output_all)
    options = utils.set_default(options, 'block_size', 500);
    options = utils.set_default(options, 'min_region_size', 9);
    options = utils.set_default(options, 'mesh_size', 100);
    options = utils.set_default(options, 'outmask_power', 3);
    options = utils.set_default(options, 'inmask_feather', 100);
    options = utils.set_default(options, 'inoutmask_ratio', 1e3);
    options = utils.set_default(options, 'prevent_identical', false);
    options = utils.set_default(options, 'glue_width', 15);
    options = utils.set_default(options, 'max_wrinke_width', 50);

    if ischar(masks)
        mask = imread(masks);
    elseif size(masks, 3) == 2
        mask = masks(:,:,1);
    else
        mask = masks;
    end
    if nargin < 4 || isempty(rid)
        rid = nan;
    end
    if nargin < 5 || isempty(outsz)
        outsz = size(mask);
    end
    if nargin < 6
        output_all = false;
    end
    if options.block_size < options.mesh_size
        error('grid spacing is smaller than the block size');
    end
    MC = struct;
    outsz = ceil(outsz / options.block_size) * options.block_size;

    bsz = options.block_size;

    Nblk_x = outsz(2) / bsz;
    Nblk_y = outsz(1) / bsz;
    [bx, by] = meshgrid(0.5:1:Nblk_x, 0.5:1:Nblk_y);
    bx = bx(:);
    by = by(:);
    Nblk = numel(bx);

    bbox = [(by - 0.5) * bsz + 1, (by + 0.5) * bsz, ...
            (bx - 0.5) * bsz + 1, (bx + 0.5) * bsz];
    block_id = [by + 0.5, bx + 0.5];
    offsets = nan(Nblk,2,'double');
    def_tiles = cell(Nblk, 1);

    % find the bbox that contains all the transformed meshes
    xys_TR = M.TR.Points;
    blk_TR = xys_TR / bsz;
    blk_TR_min = min(max(0, floor(min(blk_TR))), [Nblk_x, Nblk_y]);
    blk_TR_max = min(max(0, ceil(max(blk_TR))), [Nblk_x, Nblk_y]);
    idx_bbox = (bx > blk_TR_min(1)) & (bx < blk_TR_max(1)) & ...
               (by > blk_TR_min(2)) & (by < blk_TR_max(2));
    dxy = blk_TR_min * bsz;

    % only render the block regions covered by mesh
    xys_p = xys_TR - dxy;
    bbox_p = bbox(idx_bbox, :) - repelem(fliplr(dxy),1,2);
    outsz_p = fliplr(blk_TR_max - blk_TR_min) * bsz;
    cLists = M.TR.ConnectivityList;

    Npt = M.pt_num0;
    Ntr = M.tri_num0;
    Nrg = M.region_num;

    imght = outsz_p(1);
    imgwd = outsz_p(2);
    Ngrid_x = ceil(imgwd / options.mesh_size);
    Ngrid_y = ceil(imght / options.mesh_size);
    tx = (0:Ngrid_x) * options.mesh_size;
    ty = (0:Ngrid_y) * options.mesh_size;
    tx = tx - mean(tx(:)) + (1+imgwd)/2;
    ty = ty - mean(ty(:)) + (1+imght)/2;
    [gx, gy] = meshgrid(tx, ty);
    g_flag = false(numel(gx), 1);  % if grid points is covered by meshes

    [XX, YY] = meshgrid(1:imgwd, 1:imght);
    XY0 = cat(3, XX, YY);

    WT = zeros(imght, imgwd, 'double');
    MP = zeros(imght, imgwd, 2, 'double');
    dilat_disk = options.glue_width;

    if Nrg > 0 && ~isempty(masks)
        if size(masks, 3) == 2
            maskw = masks(:,:,2);
        else
            maskw = imfill(mask>0, 'holes') & ~ (mask>0);
        end
        Dw = bwdist(maskw);

        for k = 1:Nrg
            if any(~isnan(rid)) && all(rid ~= k)
                continue
            end
            maskg = (mask == k);
            if sum(maskg(:)) < options.min_region_size
                continue
            end
            idx1 = (k - 1) * Npt + 1;
            idx2 = k * Npt;
            xy = xys_p(idx1:idx2, :);
            idx1 = (k - 1) * Ntr + 1;
            idx2 = k * Ntr;
            cList = cLists(idx1:idx2, :) - (k-1)*Npt;
            TR1 = triangulation(double(cList), double(xy));
            [gx1, gy1, inside_grid] = local_mesh2mov_single_region(M.TR0, TR1, gx, gy);
            g_flag = g_flag | inside_grid;
            mp_x = interp2(gx, gy, gx1, (1:imgwd)', 1:imght, 'makima');
            mp_y = interp2(gx, gy, gy1, (1:imgwd)', 1:imght, 'makima');
            mp = cat(3, mp_x, mp_y);
            maskb = imdilate(maskg, ones(3)) & (mask>0) & (~maskg) & (~maskw);
            Db = bwdist(maskb);
            maskb = (Db < Dw) & (Db < dilat_disk);
            def_g = mp - XY0;
            maskd = imwarp(single(maskg), def_g) > 0.5;
            maskb = imwarp(single(maskb), def_g) > 0.5;
            maskd = max(10*maskb, maskd);
            wt = local_mask_to_wt(maskd, options.outmask_power, options.inmask_feather, ...
                                    options.inoutmask_ratio, options.max_wrinke_width);
            MP = MP + mp .* wt;
            WT = WT + wt;
        end
        MP = MP./WT;
    else
        xy = xys_p(1:Npt, :);
        cList = cLists(1:Ntr, :);
        TR1 = triangulation(double(cList), double(xy));
        [gx1, gy1, inside_grid] = local_mesh2mov_single_region(M.TR0, TR1, gx, gy);
        g_flag = g_flag | inside_grid;
        mp_x = interp2(gx, gy, gx1, (1:imgwd)', 1:imght, 'makima');
        mp_y = interp2(gx, gy, gy1, (1:imgwd)', 1:imght, 'makima');
        MP = cat(3, mp_x, mp_y);
        MP = imgaussfilt(MP, 3);
    end

    % devide the deformation field into blocks
    Nblk_p = sum(idx_bbox(:));
    offsets_p = nan(Nblk_p,2,'double');
    def_tiles_p = cell(Nblk_p,1);
    for k = 1:Nblk_p
        bx = bbox_p(k,:);
        idx0 = (gx(:) >= bx(3)-1) & (gx(:) <= bx(4)) & ...
               (gy(:) >= bx(1)-1) & (gy(:) <= bx(2));
        if all(~g_flag(idx0)) && ~output_all
            continue;
        end
        tile = MP(bx(1):bx(2), bx(3):bx(4), :);
        def_ave = nanmean(nanmean(tile,1),2);
        offsets_p(k,:) = squeeze(def_ave);
        tile = tile - def_ave;
        def_tiles_p{k} = single(tile);
    end 
    offsets(idx_bbox, :) = offsets_p;
    def_tiles(idx_bbox) = def_tiles_p;

    MC.outsz = outsz;  % [ysize, xsize]
    MC.blocksz = bsz;
    MC.Nblk = Nblk;
    MC.block_id = block_id; % [row_id, col_id]
    MC.bbox = bbox;  % [ymin, ymax, xmin, xmax]
    MC.offsets = offsets; % [x, y]
    MC.def_tiles = def_tiles;
end

function [gx1, gy1, inside_grid] = local_mesh2mov_single_region(TR0, TR1, gx, gy)
    [B, ID, oob] = elastic_mesh.cart2bary(TR1, [gy(:), gx(:)]);
    gxy1 = barycentricToCartesian(TR0,ID,B);
    gx1 = gxy1(:,1);
    gy1 = gxy1(:,2);
    gx1 = reshape(gx1, size(gx));
    gy1 = reshape(gy1, size(gy));
    inside_grid = ~oob;
end

function wt = local_mask_to_wt(mask, outmask_power, inmask_feather, inoutmask_ratio, max_wrinke_width)
    if all(~mask(:))
        wt = zeros(size(mask));
        return
    end
    [D1, idx] = bwdist(mask > 0.5);
    maskd = mask(idx);
    mask = mask > 0.5;
    D1 = max(D1, 1);
    D1(D1>max_wrinke_width) = D1(D1>max_wrinke_width) * inoutmask_ratio;
    D1 = D1.* maskd;
    D1 = 1 ./ D1 .^ outmask_power;
    D2 = bwdist(~mask);
    D2 = D2 / inmask_feather;
    % D2 = imgaussfilt(D2, 3);
    D2 = min(D2, 1) + 1000 / inoutmask_ratio;
    wt = (D2 .* mask) + (D1 .* (1 - mask)) / inoutmask_ratio;
end