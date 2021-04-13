function maskr = partition_section_extend_wrinkle_unique(maskw, masko, options)
    if nargin < 3
        options = [];
    end
    options = utils.set_default(options, 'multiplier', 2);
    options = utils.set_default(options, 'minhalflen', 2);
    options = utils.set_default(options, 'dyxcutoff', 0.2);
    options = utils.set_default(options, 'angle_thresh', 35);
    LBL_O = 3;
    LBL_W = 2;
    LBL_A = 1;
    if all(~maskw(:))
        maskr = ones(size(maskw), 'uint8');
        maskr(masko) = 0;
    end
    [endpt, L] = mask_gen.detect_end_pt(maskw, options);
    mask = LBL_O*masko + LBL_W*maskw;
    Npt = size(endpt, 1);
    % pt_flag = false(Npt, 1);
    pps = nchoosek(1:Npt, 2);
    to_connect = false(size(pps,1),1);
    B = zeros(size(pps,1),1);
    for k = 1:size(pps,1)
        idx1 = pps(k, 1);
        idx2 = pps(k, 2);
        if endpt.id(idx1) == endpt.id(idx2)
            continue;
        end
        vyx = diff(endpt.yx([idx1, idx2],:),1,1);
        vyx = vyx / sum(vyx(:).^2).^0.5;
        b1 = sum(-vyx.*endpt.dyx(idx1,:));
        b2 = sum(vyx.*endpt.dyx(idx2,:));
        B(k) = min(b1, b2);
        if b1 > cosd(options.angle_thresh) && b2 > cosd(options.angle_thresh)
            to_connect(k) = true;
            % pt_flag([idx1, idx2]) = true;
        end
    end
    pps = pps(to_connect,:);
    B = B(to_connect);
    [~, idx] = sort(B,'descend');
    pps = pps(idx, :);
    % ensure one way
    to_connect = true(size(pps,1),1);
    for k = 1:(size(pps,1)-1)
        if ~to_connect(k)
            continue
        end
        pp = pps(k,:);
        t1 = ~any(pps == pp(1) | pps == pp(2), 2);
        to_connect(k+1:end,:) = to_connect(k+1:end,:) & t1(k+1:end,:);
    end
    pps = pps(to_connect,:);

    start_yx = endpt.yx(pps(:,1),:);
    end_yx =  endpt.yx(pps(:,2),:);
    id1 = endpt.id(pps(:,1));
    id2 = endpt.id(pps(:,2));
    indx = mask_gen.add_lines_enpts(size(maskw), start_yx, end_yx, true);
    Lindx = L(indx);
    crossed_flag = all((Lindx == 0) | (Lindx == id1) | (Lindx == id2),2);
    indx = indx(crossed_flag,:);
    pps = pps(crossed_flag, :);
    % prevent crossing
    start_yx = endpt.yx(pps(:,1),:);
    end_yx =  endpt.yx(pps(:,2),:);
    kept_flag = true(size(indx,1),1);
    for k = 2:size(indx,1)
        for k1 = 1:(k-1)
            if ~kept_flag(k1)
                continue;
            end
            if local_if_cross(start_yx([k,k1],:),end_yx([k,k1],:))
                kept_flag(k) = false;
                break
            end
        end
    end
    indx = indx(kept_flag,:);
    pps = pps(kept_flag,:);
    
    pt_flag = ismember((1:Npt)', pps(:));
    maskt = zeros(size(maskw));
    maskt(indx) = LBL_A;
    mask = max(mask, maskt);
    if all(pt_flag)
        return
    end
    indx2 = mask_gen.add_line_direct(size(maskw), endpt.yx(~pt_flag,:), -endpt.dyx(~pt_flag,:), false);
    id = endpt.id(~pt_flag);
    line_num = cellfun(@numel, indx2);
    [~, idxt] = sort(line_num);
    indx2 = indx2(idxt);
    id = id(idxt);
    for k = 1:numel(indx2)
        maskt = zeros(size(maskw));
        maskt(indx2{k}) = LBL_A;
        maskn = (mask > 0) & (mask < LBL_O) & (L ~= id(k));
        if any(maskt(:) & maskn(:))
            maskt = LBL_A * imreconstruct(L==id(k), maskt & ~maskn);
        end
        mask = max(mask, maskt); 
    end
    mask = max(mask, LBL_A*imdilate(mask>0,ones(3)));
    maskr = bwlabel(mask==0, 4);
    [~,idx] = bwdist(mask==0);
    maskr = uint8(maskr(idx) .* (mask < LBL_W));
end

function c = local_if_cross(sttyx, endyx)
    AC = sttyx(2,:) - sttyx(1,:);
    BD = endyx(2,:) - endyx(1,:);
    BC = sttyx(2,:) - endyx(1,:);
    AD = endyx(2,:) - sttyx(1,:);
    c1 = local_cross_product(AC, AD) * local_cross_product(BC, BD) < 0;
    c2 = local_cross_product(AC, BC) * local_cross_product(AD, BD) < 0;
    c = c1 & c2;
end

function p = local_cross_product(x, y)
    p = x(:,1) .* y(:,2) - x(:,2) * y(:,1);
end