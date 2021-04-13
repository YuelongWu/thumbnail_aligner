function blk = retrieve_same_block(img, kps0, kd, rid)
    if nargin < 4
        rid = [];
    end
    if isfield(kps0, 'block_id')
        nb_size = kps0.nb_size(kd);
        blk = PMCC.divide_img_to_block(img, nb_size, 0);
        block_id = kps0.block_id{kd};
        if ~isempty(rid)
            region_id = kps0.region_id{kd};
            idxt = ismember(region_id, rid);
            block_id = block_id(idxt); 
        end
        blk = blk(:,:,block_id);
    else
        imgwd = size(img,2);
        imght = size(img,1);
        yx = kps0.yx{kd};
        nb_size = kps0.nb_size(kd);
        if ~isempty(rid)
            region_id = kps0.region_id{kd};
            idxt = ismember(region_id, rid);
            yx = yx(idxt, :);
        end
        yy = permute(yx(:,1), [2,3,1]);
        xx = permute(yx(:,2), [2,3,1]);
        wsz = round(nb_size / 2);
        [xg, yg] = meshgrid(-wsz:1:wsz, -wsz:1:wsz);
        xx = xx + xg;
        yy = yy + yg;
        xx = min(max(xx, 1),imgwd);
        yy = min(max(yy, 1),imght);
        indx = sub2ind([imght, imgwd], yy, xx);
        blk = img(indx);
    end
end
