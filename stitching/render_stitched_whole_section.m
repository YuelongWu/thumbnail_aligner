function [IMGt, outmask, Ms] = render_stitched_whole_section(imglist, Ms, options)
    options = utils.set_default(options, 'rotate_image', false);
    options = utils.set_default(options, 'block_sz', 0);
    options = utils.set_default(options, 'mesh_sz', 100);
    options = utils.set_default(options, 'blend_range', 500);
    options = utils.set_default(options, 'weight_power', 1);
    options = utils.set_default(options, 'blend_method', 'mean');
    options = utils.set_default(options, 'interp', 'linear1');
    options = utils.set_default(options, 'fillval', 0);
    Nimg = numel(imglist);

    if options.rotate_image
        thetas = zeros(Nimg, 1);
        for t = 1:Nimg
            M = Ms{t};
            [~, R] = geometries.fit_affine(M.TR0.Points, M.TR.Points);
            thetas(t) = atan2(R(1,2), R(1,1));
        end
        theta = mean(thetas);
        R = [cos(theta), sin(theta); -sin(theta), cos(theta)];
    end

    xy_min = nan(Nimg,2);
    xy_max = nan(Nimg,2);
    for t = 1:Nimg
        M = Ms{t};
        if options.rotate_image
            M.TR.Points = M.TR.Points * R;
            Ms{t} = M;
        end
        xy_min(t,:) = min(M.TR.Points, [], 1);
        xy_max(t,:) = max(M.TR.Points, [], 1);
    end
    offst = min(xy_min);
    outsz = fliplr(ceil(max(xy_max) - offst));
    if options.block_sz > 0
        outsz0 = ceil(outsz / block_sz) * block_sz;
        offst0 = offst - fliplr(outsz0 - outsz)/2;
        outsz = outsz0;
        offst = offst0;
    end
    IMGt = zeros(outsz);
    WT = zeros(outsz);
    crnt_sz = [nan, nan];
    for k = 1:numel(imglist)
        % if isfield(imglist(k), 'img') && ~isempty(imglist(k).img)
        %     tile = single(imglist(k).tile);
        % else
        %     tile = single(imread([imglist(k).folder, filesep, imglist(k).name]));
        % end
        tile = single(imread([imglist(k).folder, filesep, imglist(k).name]));
        M = Ms{k};
        M.TR.Points = M.TR.Points - offst;
        Ms{k} = M;
        tilesz = [size(tile,1), size(tile,2)];
        if any(tilesz ~= crnt_sz)
            wt = local_generate_mask(tilesz, options);
            crnt_sz = tilesz;
        end
        [D, bbox] = local_tile_deformation(M, options.mesh_sz, [1,outsz(1),1,outsz(2)]);
        tile_n = imwarp(tile, D, options.interp);
        wt_n = imwarp(wt, D, options.interp);
        if strcmpi(options.blend_method, 'mean')
            IMGt(bbox(1):bbox(2),bbox(3):bbox(4)) = ...
                IMGt(bbox(1):bbox(2),bbox(3):bbox(4)) + tile_n .* wt_n;
            WT(bbox(1):bbox(2),bbox(3):bbox(4)) = ...
                WT(bbox(1):bbox(2),bbox(3):bbox(4)) + wt_n;
        else
            wt0 = WT(bbox(1):bbox(2),bbox(3):bbox(4));
            idxt = wt0(:) <= wt_n(:);
            imgt = IMGt(bbox(1):bbox(2),bbox(3):bbox(4));
            IMGt(bbox(1):bbox(2),bbox(3):bbox(4)) = ...
                imgt(:) .* (1-idxt(:)) +  tile_n(:) .* wt_n(:) .* idxt(:);
            WT(bbox(1):bbox(2),bbox(3):bbox(4)) = max(wt_n, wt0);
        end
    end
    IMGt = IMGt ./ WT;
    outmask = isnan(IMGt);
    % if options.fillval == 0
    %     replcval = 1;
    % else
    %     replcval = options.fillval - 1;
    % end
    % IMGt((~outmask) & (IMGt == options.fillval)) = replcval;
    IMGt(outmask) = options.fillval;
    IMGt = uint8(IMGt);
end

function wt = local_generate_mask(imgsz, options)
    if options.weight_power == 0
        wt = ones(imgsz);
        return
    end
    wt = zeros(imgsz);
    wt([1,end],:) = 1;
    wt(:,[1,end]) = 1;
    wt = bwdist(wt) + 1;
    wt = wt / options.blend_range;
    wt = min(wt,1);
    wt = wt .^ options.weight_power;
end

function [D, bbox] = local_tile_deformation(M, meshsz, yx_limit)
    % yx_limit, bbox: [ymin, ymax, xmin, xmax]
    if nargin < 3
        yx_limit = [];
    end
    XY1 = M.TR.Points;
    xymin = floor(min(XY1)) + 1;
    xymax = ceil(max(XY1));
    xymin = max(xymin, [yx_limit(3), yx_limit(1)]);
    xymax = min(xymax, [yx_limit(4), yx_limit(2)]);
    xy = XY1 - xymin + 1;
    imgsz = fliplr(xymax - xymin + 1);
    cList = M.TR.ConnectivityList;
    TR1 = triangulation(double(cList), double(xy));
    MC = elastic_mesh.mesh_to_movingcoord(M.TR0, TR1, imgsz, meshsz);
    [gxx, gyy] = meshgrid(1:imgsz(2), 1:imgsz(1));
    D = MC - cat(3, gxx, gyy);
    bbox = [xymin(2), xymax(2), xymin(1), xymax(1)];
end

