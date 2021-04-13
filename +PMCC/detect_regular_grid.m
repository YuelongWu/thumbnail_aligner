function kps = detect_regular_grid(img0, options, mask)
    kps = struct;
    options = utils.set_default(options, 'min_dis', 100);
    options = utils.set_default(options, 'nb_size', 201);
    options = utils.set_default(options, 'fft', false);
    options = utils.set_default(options, 'pad_des', true);

    if (nargin < 3) || isempty(mask)
        mask = 1;
    end
    if ischar(mask)
        mask = imread(mask);
    end
    mask = cast(mask, class(img0));
    rgn_id = unique(mask(mask > 0));
    N_rgn = numel(rgn_id);

    min_dis = options.min_dis;
    nb_size = options.nb_size;
    [~, ia_g] = unique([min_dis(:), nb_size(:)], 'rows');
    min_dis = min_dis(ia_g);
    nb_size = nb_size(ia_g);

    kps.min_dis = min_dis;
    kps.nb_size = nb_size;

    N_dis = numel(min_dis);
    YXs = cell(N_rgn, N_dis);
    VALs = cell(N_rgn, N_dis);
    RIDs = cell(N_rgn, N_dis);
    DESs = cell(N_rgn, N_dis);

    imght = size(img0,1);
    imgwd = size(img0,2);
    imgsz = [imght, imgwd];
    [xx0, yy0] = meshgrid(1:imgwd, 1:imght);
    for kr = 1:N_rgn
        maskt = mask == rgn_id(kr);
        img = img0 .* maskt;
        if all(img(:) <= 0)
            continue
        end
        cntr_x = mean(xx0(maskt));
        cntr_y = mean(yy0(maskt));
        cntr = [cntr_y, cntr_x];
        for kd = 1:N_dis
            yx_grid = local_gen_eqtriang_mesh(maskt, min_dis(kd),0, cntr);
            indx = sub2ind(imgsz, yx_grid(:,1), yx_grid(:,2));
            idx = mask(indx) == rgn_id(kr);
            if sum(idx(:)) < 3
                D = bwdist(mask == rgn_id(kr));
                Didx = D(indx);
                [~, idx] = mink(Didx, 3);
            end
            yx = yx_grid(idx, :);
            if isempty(yx)
                continue
            end
            indx = indx(idx);
            val = img(indx);
            [val, idx] = sort(val, 'descend', 'MissingPlacement', 'last','ComparisonMethod','abs');
            yx = yx(idx, :);

            yy = permute(yx(:,1), [2,3,1]);
            xx = permute(yx(:,2), [2,3,1]);
            wsz = round(nb_size(kd) / 2);
            [xg, yg] = meshgrid(-wsz:1:wsz, -wsz:1:wsz);
            xx = xx + xg;
            yy = yy + yg;
            xx = min(max(xx, 1),imgwd);
            yy = min(max(yy, 1),imght);
            indx = sub2ind([imght, imgwd], yy, xx);
            des = img(indx);
            des1d = reshape(des, size(des,1)*size(des,2), size(des,3));
            valid_idx = std(des1d, 1) > 0;
            if options.pad_des == 1
                des = padarray(des, floor([size(des,1), size(des,2), 0]/2), 'both');
            elseif options.pad_des == 2
                pdwd = floor([size(des,1), size(des,2), 0]/2);
                if pdwd(1) > 2 * pdwd(2)
                    pdwd(1) = 0;
                elseif 2 * pdwd(1) < pdwd(2)
                    pdwd(2) = 0;
                end
                des = padarray(des, pdwd, 'both');
            end
            YXs{kr, kd} = yx(valid_idx,:);
            VALs{kr, kd} = val(valid_idx);
            RIDs{kr, kd} = rgn_id(kr) * ones(numel(val(valid_idx)), 1);
            DESs{kr, kd} = des(:,:,valid_idx);
        end
    end
    kps.yx = cell(N_dis,1);
    kps.val = cell(N_dis,1);
    kps.region_id = cell(N_dis,1);
    kps.des = cell(N_dis,1);
    if numel(options.fft) == 1
        kps.ffted = repmat(options.fft, N_dis, 1);
    else
        kps.ffted = options.fft(ia_g);
    end
    for kd = 1:N_dis
        kps.yx{kd} = vertcat(YXs{:, kd});
        kps.val{kd} = vertcat(VALs{:, kd});
        kps.region_id{kd} = vertcat(RIDs{:,kd});
        if kps.ffted(kd)
            kps.des{kd} = fft2(cat(3, DESs{:, kd}));
        else
            kps.des{kd} = cat(3, DESs{:, kd});
        end
    end
end

function yx = local_gen_eqtriang_mesh(mask, sideL, extend, cntr)
    if nargin < 4
        cntr = [];
    end
    if nargin < 3 || isempty(extend)
        extend = 0;
    end
    theta = 0;
    imgsz = size(mask);
    imght = size(mask,1);
    imgwd = size(mask,2); 
    vec1 = [sind(theta), cosd(theta)] * sideL;
    vec2 = [sind(theta + 60), cosd(theta + 60)] * sideL;
    vec1n = vec1 - vec2 * dot(vec1, vec2) / dot(vec2, vec2);
    vec2n = vec2 - vec1 * dot(vec1, vec2) / dot(vec1, vec1);
    sideLn2 = dot(vec1n, vec1n);
    gridht = imght;
    gridwd = imgwd;
    if extend > 0
        gridht = gridht + 2 * sideL * (extend + 0.25);
        gridwd = gridwd + 2 * sideL * (extend + 0.25);
    end
    corners = [0, 0; gridht, 0; 0, gridwd; gridht, gridwd];
    cornercoord = corners * [vec1n(:), vec2n(:)] / sideLn2;
    coordrange = ceil(range(cornercoord, 1));
    [vv1, vv2] = meshgrid(1:coordrange(1), 1:coordrange(2));
    V = [vv1(:),vv2(:)];

    pts = V * [vec1; vec2];
    pts = pts - mean(pts, 1) + [gridht, gridwd]/2;
    if ~isempty(cntr)
        tmpdis = sum((pts - cntr).^2, 2);
        [~, tmpidx] = min(tmpdis);
        pts = pts - pts(tmpidx, :) + cntr;
    end
    sbs = min(max(1, round(pts)), imgsz);
    ind = sub2ind(imgsz, sbs(:,1), sbs(:,2));
    idx0 = (pts(:, 1) <= imght) &  (pts(:, 1) >= 0) & ...
        (pts(:, 2) <= imgwd) &  (pts(:, 2) >= 0);
    idx1 = mask(ind) > 0;
    idx = idx0 & idx1;
    pts0 = pts(idx, :);
    if extend > 0
        pts_outside = pts(~idx, :);
        dis = pdist2(pts0, pts_outside, 'euclidean', 'Smallest', 1);
        idxt = dis < sideL * (extend + 0.1);
        yx = [pts0; pts_outside(idxt,:)];
    else
        yx = pts0;
    end
    yx = round(yx);
    yx = max(1, min(yx, imgsz));
end

