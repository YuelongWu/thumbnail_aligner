function [D, WT] = subregion_movingcoord(M, masks, options, rid, outsz)
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
    if nargin < 5
        outsz = size(mask);
    end
    mask = double(mask);
    if nargout > 1
        maskexp = mask_gen.expand_mask(mask);
    end
    xys = M.TR.Points;
    cLists = M.TR.ConnectivityList;
    Npt = M.pt_num0;
    Ntr = M.tri_num0;
    Nrg = M.region_num;
    imght = outsz(1);
    imgwd = outsz(2);
    % xys(:,1) = xys(:,1) + 1000;
    % imgwd = imgwd + 2000;
    WTs = zeros(imght,imgwd,1,Nrg);
    Ds = zeros(imght,imgwd,2,Nrg);
    [XX, YY] = meshgrid(1:imgwd, 1:imght);
    XY0 = cat(3, XX, YY);
    dilat_disk = options.glue_width;
    MASK = ones(imght,imgwd,1,Nrg);
    if size(masks,3) == 2
        maskw = masks(:,:,2);
    else
        maskw = imfill(mask>0, 'holes') & ~ (mask>0);
    end
    Dw = bwdist(maskw);
    % maskw = imdilate(maskw, strel('disk',round(dilat_disk)));
    % maskp = mask > 0;
    for k = 1:Nrg
        if all(~isnan(rid)) && all(rid ~= k)
             continue
        end
        idx1 = (k - 1) * Npt + 1;
        idx2 = k * Npt;
        xy = xys(idx1:idx2,:);
        if options.prevent_identical && all(xy(:) == M.TR0.Points(:))
            continue
        end
        idx1 = (k - 1) * Ntr + 1;
        idx2 = k * Ntr;
        cList = cLists(idx1:idx2, :) - (k-1)*Npt;
        TR1 = triangulation(double(cList), double(xy));
        XY1 = elastic_mesh.mesh_to_movingcoord(M.TR0, TR1, outsz, options.mesh_size);
        Ds(:,:,:,k) = XY1;
        Dt = XY1 - XY0;
        maskt = mask == k;
        maskb = imdilate(maskt,ones(3)) & (mask > 0) & (~maskt) & (~maskw);
        % maskb = maskb & imerode(maskp, strel('disk', round(dilat_disk/2)));
        Db = bwdist(maskb);
        maskb = (Db < Dw) & (Db < dilat_disk);
        % maskb = imdilate(maskb, strel('disk', dilat_disk));

        maskd = imwarp(single(maskt), Dt, 'linear') > 0.5;
        maskb = imwarp(single(maskb), Dt) > 0.5;
        maskd = max(10*maskb, maskd);
        if nargout > 1
            maskdexp = imwarp(single((maskt > 0.5) | (maskexp == k)), Dt, 'linear') > 0.5;
            MASK(:,:,k) = maskdexp;
        end
        wt = mask_to_wt(maskd, options.outmask_power, options.inmask_feather, options.inoutmask_ratio, options.max_wrinke_width);
        WTs(:,:,:,k) = wt;
    end
    D = sum(Ds .* WTs, 4) ./ sum(WTs, 4);
    WT = sum(MASK, 4);
end

function wt = mask_to_wt(mask, outmask_power, inmask_feather, inoutmask_ratio, max_wrinke_width)
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
    D2 = min(D2, 1) + 10 / inoutmask_ratio;
    wt = (D2 .* mask) + (D1 .* (1 - mask)) / inoutmask_ratio;
end
