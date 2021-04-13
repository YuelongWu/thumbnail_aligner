function kps = filter_kps(kps0, kd, rid)
    if nargin < 3
        rid = [];
    end
    kps = struct;
    kps.nb_size = kps0.nb_size(kd);
    yx = kps0.yx{kd};
    val = kps0.val{kd};
    if isfield(kps0,'des')
        des = kps0.des{kd};
    end
    region_id = kps0.region_id{kd};
    if ~isempty(rid)
        idxt = ismember(region_id, rid);
        yx = yx(idxt,:);
        val = val(idxt,:);
        if isfield(kps0,'des')
            des = des(:,:,idxt);
        end
        region_id = region_id(idxt);
    else
        idxt = true(size(kps0.region_id{kd}));
    end
    kps.yx = yx;
    kps.val = val;
    if isfield(kps0,'des')
        kps.des = des;
    end
    kps.region_id = region_id;
    kps.ffted = kps0.ffted(kd);
    if isfield(kps0, 'block_id')
        block_id = kps0.block_id{kd};
        block_id = block_id(idxt);
        kps.block_id = block_id;
    end
    if isfield(kps0, 'min_dis')
        kps.min_dis = kps0.min_dis(kd);
    end
end
