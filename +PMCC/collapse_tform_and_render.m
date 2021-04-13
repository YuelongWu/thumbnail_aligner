function [TFORMS, IMGt] = collapse_tform_and_render(IMG, TFORMS, gridsz, maskr, rids)
    if nargin < 5 || isempty(rids)
        rids = unique(maskr(maskr>0));
    end
    imgtp = class(IMG);
    if ~strcmpi(imgtp, 'single')
        IMG = single(IMG);
    end
    methd = 'linear';
    imght = size(maskr, 1);
    imgwd = size(maskr, 2);
    Nrg = numel(rids);
    Ntf = numel(TFORMS);
    bboxes = zeros(Nrg, 4);
    PTPAIRS = cell(Nrg, 1);

    if ~isempty(IMG)
        xg1 = 1 : gridsz : imgwd;
        xg2 = xg1 - 1 + gridsz;
        xg = (xg1 + xg2) / 2;
        yg1 = 1 : gridsz : imght;
        yg2 = yg1 - 1 + gridsz;
        yg = (yg1 + yg2) / 2;


        imght0 = ceil(imght / gridsz) * gridsz;
        imgwd0 = ceil(imgwd / gridsz) * gridsz;
        IMGt = zeros(imght0, imgwd0, class(IMG));
        maskrt = zeros(imght0, imgwd0, class(maskr));
        maskrt(1:imght, 1:imgwd) = maskr;
        maskr = maskrt;
    else
        IMGt = [];
    end
    rid2_defined = false;
    for kr = 1:Nrg
        rid = rids(kr);
        % crop mask to save some time
        mask = maskr == rid;

        if ~isempty(IMG)
            [~, bbox] = mask_gen.crop_mask(mask, mask);
            bboxes(kr, :) = bbox;
            id_y1 = find(yg1 <= bbox(1), 1, 'last');
            id_y2 = find(yg2 >= bbox(3), 1, 'first');
            id_x1 = find(xg1 <= bbox(2), 1, 'last');
            id_x2 = find(xg2 >= bbox(4), 1, 'first');
            xgr = xg(id_x1:id_x2);
            ygr = yg(id_y1:id_y2);
            [xx0, yy0] = meshgrid(xgr, ygr);
            yx1 = [yy0(:), xx0(:)];
            gsh = size(xx0);
        else
            yx1 = gridsz;
            yx1_rnd = round(yx1);
            idxt = all(yx1_rnd > 0 & yx1_rnd <= [imght, imgwd],2);
            yx1 = yx1(idxt,:);
            yx1_rnd = yx1_rnd(idxt,:);
            ind = sub2ind(size(mask), yx1_rnd(:,1), yx1_rnd(:,2));
            idxt = mask(ind);
            yx1 = yx1(idxt,:);
        end
        yx2 = yx1;
        conf_defined = false;
        conf = ones(size(yx2,1),1, 'single');
        rotation = zeros(size(yx2,1), 1, 'single');
        for kt = Ntf:-1:1
            tfm = TFORMS{kt};
            if isempty(tfm)
                if ~rid2_defined
                    rid2 = 1;
                end
                continue
            end
            if isstruct(tfm) % affine
                idxt = find(tfm.rid(:,1) == rid, 1);
                if ~rid2_defined
                    rid2 = tfm.rid(1,2);
                    rid2_defined = true;
                end
                if ~isempty(idxt)
                    A0 = tfm.A(:,:,idxt);
                    yx2 = yx2 * A0(1:2,1:2) + A0(3,1:2);
                end
            end
            if istable(tfm)
                idxt = tfm.region_id(:, 1) == rid;
                if ~rid2_defined
                    rid2 = tfm.region_id(1,2);
                    rid2_defined = true;
                end
                if any(idxt)
                    ptpairs = tfm(idxt, :);
                    yx1t = double(ptpairs.yx1);
                    yx2t = double(ptpairs.yx2);
                    rnk = rank(geometries.padones(yx2t));
                    if rnk <= 2 || size(yx1t,1) <= 3
                        FA = geometries.fit_affine(yx2t, yx1t);
                        yx2 = yx2 * FA(1:2,1:2) + FA(3,1:2);
                    else
                        Fx = scatteredInterpolant(yx1t, yx2t(:,2), methd);
                        Fy = scatteredInterpolant(yx1t, yx2t(:,1), methd);
                        yx2_buf = [Fy(yx2), Fx(yx2)];
                        if isempty(yx2_buf)
                            FA = geometries.fit_affine(yx2t, yx1t);
                            yx2 = yx2 * FA(1:2,1:2) + FA(3,1:2);
                        else
                            yx2 = yx2_buf;
                        end
                    end
                    if ~conf_defined
                        Fc = scatteredInterpolant(yx1t, double(ptpairs.conf), 'nearest', 'nearest');
                        conf = Fc(yx2);
                        if isempty(conf)
                            conf = nanmean(ptpairs.conf) * ones(size(yx2,1), 1, 'single');
                        end
                        conf_defined = true;
                    end
                end
            end
        end
        region_id = repmat(single([rid, rid2]), size(yx1,1),1);
        PTPAIRS{kr} = table(yx1, yx2, region_id, conf, rotation);
        if ~isempty(IMG)
            XM = reshape(yx2(:,2), gsh);
            YM = reshape(yx2(:,1), gsh);
            [XF, YF] = meshgrid(xg(1:size(XM,2)), yg(1:size(XM,1)));
            XD = XM - XF;
            YD = YM - YF;
            if any(isnan(XD(:)))
                maskX = isnan(XD);
                [~, idxX] = bwdist(~maskX);
                XD = XD(idxX);
            end
            if any(isnan(YD(:)))
                maskY = isnan(YD);
                [~, idxY] = bwdist(~maskY);
                YD = YD(idxY);
            end
            Dd = cat(3, XD, YD);
            D = imresize(Dd, gridsz,'bicubic');
            IMGr = imwarp(IMG, D, 'linear','FillValues', 0);
            maskt = mask(yg1(id_y1):yg2(id_y2), xg1(id_x1):xg2(id_x2));
            IMGt(yg1(id_y1):yg2(id_y2), xg1(id_x1):xg2(id_x2)) = ...
                IMGt(yg1(id_y1):yg2(id_y2), xg1(id_x1):xg2(id_x2)) + IMGr .* maskt;
        end
    end
    idxt = cellfun(@isempty, PTPAIRS);
    PTPAIRS = PTPAIRS(~idxt);
    TFORMS = vertcat(PTPAIRS{:});
    if ~isempty(IMG)
        IMGt = IMGt(1:imght, 1:imgwd);
        if ~strcmpi(imgtp, 'single')
            IMGt = cast(IMGt, imgtp);
        end
    end
end
