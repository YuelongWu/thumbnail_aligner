function ptpairs = block_match(kps1, kps2, options, ptpairs0)
    % match smaller neighbor in kps1 to larger neighbor in kps2
    if nargin > 3 && ~isempty(ptpairs0)
        [rgpairs, Apm, Rpm] = geometries.prematch2affine(ptpairs0);
        Np = size(rgpairs, 1);
    else
        rid1 = unique(kps1.region_id{1});
        rid2 = unique(kps2.region_id{2});
        rgpairs = combvec(rid1(:)', rid2(:)')';
        Np = size(rgpairs, 1);
        Apm = repmat(eye(3),1,1,Np);
        Rpm = repmat(eye(3),1,1,Np); 
    end
    
    options = utils.set_default(options, 'rot_thresh', 5);
    options = utils.set_default(options, 'search_range', 300);
    options = utils.set_default(options, 'mask_disk', 5);
    options = utils.set_default(options, 'min_conf', 0.5);

    yx1s = kps1.yx{1};
    rid1s = kps1.region_id{1};
    yx2s = kps2.yx{2};
    rid2s = kps2.region_id{2};

    YX1 = cell(Np, 1);
    YX2 = cell(Np, 1);
    RID = cell(Np, 1);
    CONF = cell(Np, 1);
    ROT = cell(Np, 1);
    for k = 1:Np
        rgpr = rgpairs(k, :);
        idx1 = find(rid1s == rgpr(1));
        idx2 = find(rid2s == rgpr(2));
        yx1t = yx1s(idx1, :);
        yx1tt = yx1s(idx1, :) * Apm(1:2,1:2,k) + Apm(3,1:2,k);
        yx2t = yx2s(idx2, :);
        dis = pdist2(yx1tt, yx2t, 'chebychev');
        inrange = dis <= options.search_range;
        if all(~inrange(:))
            continue
        end
        [idx1t, idx2t] = find(inrange);
        idx1tt = idx1(idx1t(:));
        idx2tt = idx2(idx2t(:));
        des1 = kps1.des{1}(:,:,idx1tt);
        des2 = kps2.des{2}(:,:,idx2tt);

        rot = rad2deg(atan2(Rpm(1,2,k), Rpm(1,1,k)));
        dyx_expected = yx1tt(idx1t,:) - yx2t(idx2t,:);
        if abs(rot) > options.rot_thresh
            des1 = imrotate(des1, rot);
        end
        [dyx, conf] = PMCC.xcorr_fft(des1, des2, ...
            options.mask_disk, dyx_expected, [kps1.ffted(1), kps2.ffted(2)], options.search_range);
        yx1 = yx1t(idx1t, :);
        yx2 = yx2t(idx2t, :) + dyx;
        confidx = conf(:) > options.min_conf;
        YX1{k} = yx1(confidx,:);
        YX2{k} = yx2(confidx,:);
        RID{k} = repmat(rgpr, sum(confidx(:)), 1);
        CONF{k} = conf(confidx);
        ROT{k} = rot * ones(sum(confidx(:)), 1);
    end

    yx1 = cat(1, YX1{:});
    yx2 = cat(1, YX2{:});
    region_id = cat(1, RID{:});
    conf = cat(1, CONF{:});
    rotation = cat(1, ROT{:});
    ptpairs = table(yx1, yx2, region_id, conf, rotation);
end
