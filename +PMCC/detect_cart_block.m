function kps = detect_cart_block(img0, options, mask, return_des)
    kps = struct;
    options = utils.set_default(options, 'nb_sz', 300);
    options = utils.set_default(options, 'fft', false);
    options = utils.set_default(options, 'pad_des', false);

    if nargin < 4
        return_des = true;
    end
    if (nargin < 3) || isempty(mask)
        mask = ones(size(img0),'single');
    end
    if ischar(mask)
        mask = imread(mask);
    end
    mask = cast(mask, class(img0));
    rgn_id = unique(mask(mask > 0));
    rgn_id = sort(rgn_id);
    rgn_int = single(rgn_id);
    rgn_int = [min(rgn_int(:))-1; rgn_int(:); max(rgn_int(:))+1];
    rgn_int = (rgn_int(2:end) + rgn_int(1:end-1))/2;
    N_rgn = numel(rgn_id);
    

    [nb_size, ia_g] = unique(options.nb_size);
    kps.nb_size = nb_size;
    N_dis = numel(nb_size);

    if numel(options.fft) == 1
        kps.ffted = repmat(options.fft, N_dis, 1);
    else
        kps.ffted = options.fft(ia_g);
    end

    kps.yx = cell(N_dis,1);
    kps.val = cell(N_dis,1);
    kps.region_id = cell(N_dis,1);
    kps.block_id = cell(N_dis,1);
    if return_des
        kps.des = cell(N_dis,1);
    end

    for kd = 1:N_dis
        bsz = nb_size(kd);
        [maskb, yx] = PMCC.divide_img_to_block(mask, bsz, 0);
        if return_des
            blk = PMCC.divide_img_to_block(img0, bsz, 0);
        end
        Nblk = size(yx, 1);
        maskb = reshape(maskb, bsz * bsz, Nblk);
        rgcnts = histc(maskb, rgn_int);
        maskb = reshape(maskb, bsz, bsz, Nblk);
        YXs = cell(N_rgn,1);
        VALs = cell(N_rgn,1);
        RIDs = cell(N_rgn,1);
        DESs = cell(N_rgn,1);
        BIDs = cell(N_rgn,1);
        for kr = 1:N_rgn
            cnt = rgcnts(kr, :);
            idxt = find(cnt > 0);
            YXs{kr} = yx(idxt, :);
            cntt = cnt(idxt);
            VALs{kr} = cntt(:);
            RIDs{kr} = single(rgn_id(kr)) * ones(numel(idxt(:)),1);
            BIDs{kr} = idxt(:);
            if return_des
                 des0 = blk(:,:,idxt) .* cast(maskb(:,:,idxt) == rgn_id(kr), class(blk));
                 if options.pad_des
                    des0 = padarray(des0, floor([size(des0,1), size(des0,2), 0]/2), 'both');
                 end
                 DESs{kr} = des0;
            end
        end
        kps.yx{kd} = vertcat(YXs{:});
        kps.val{kd} = vertcat(VALs{:});
        kps.region_id{kd} = vertcat(RIDs{:});
        kps.block_id{kd} = vertcat(BIDs{:});
        if return_des
            if kps.ffted(kd)
                kps.des{kd} = fft2(cat(3,DESs{:}));
            else
                kps.des{kd} = cat(3, DESs{:});
            end
        end
    end
end