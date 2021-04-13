function [endpt, L] = detect_end_pt(mask, options)
    if nargin < 2
        options = [];
    end
    options = utils.set_default(options, 'multiplier', 3);
    options = utils.set_default(options, 'minhalflen', 2);
    options = utils.set_default(options, 'dyxcutoff', 1);
    CC = bwconncomp(mask);
    L = labelmatrix(CC);
    imgsz = size(mask);
    maskctr = false(imgsz);
    De = bwdist(~mask);
    wd = zeros(CC.NumObjects, 1);
    for k = 1:CC.NumObjects
        [yy, xx] = ind2sub(imgsz, CC.PixelIdxList{k});
        dis2  = (xx - mean(xx)).^2 + (yy - mean(yy(:))).^2;
        [~, idx] = min(dis2);
        maskctr(yy(idx), xx(idx)) = true;
        wd(k) = 4 * mean(De(CC.PixelIdxList{k}));
    end
    Dg = bwdistgeodesic(mask, maskctr);
    
    yx = cell(CC.NumObjects,1);
    id = cell(CC.NumObjects,1);
    dyx = cell(CC.NumObjects,1);
    for k = 1:CC.NumObjects
        idxt = CC.PixelIdxList{k};
        [yy, xx] = ind2sub(imgsz, idxt);
        Dgt = nan(imgsz);
        Dgt(idxt) = Dg(idxt);
        seL = max(3, ceil(options.multiplier * wd(k)));
        Dgtd = imdilate(imdilate(max(Dgt, 0), ones(seL,1)), ones(1,seL));
        maskep = (Dgtd == Dgt) & (Dgt > options.minhalflen * wd(k));
        % seLt = max(3, ceil(wd(k)));
        % maskep = imdilate(maskep, ones(seLt)) & (Dgt>0);
        CCt = bwconncomp(maskep);
        yxt = nan(CCt.NumObjects, 2);
        dyxt = nan(CCt.NumObjects, 2);
        tokeep = true(CCt.NumObjects,1);
        for t = 1:CCt.NumObjects
            [yyt, xxt] = ind2sub(imgsz, CCt.PixelIdxList{t});
            yyt = mean(yyt(:));
            xxt = mean(xxt(:));
            if any(pdist2([yyt,xxt], yxt) <  2 * wd(k))
                tokeep(t) = false;
                continue
            end
            yxt(t,:) = round([yyt, xxt]);
            
            wt = ((yy(:) - yyt).^2 + (xx(:) - xxt).^2).^0.5;
            wt = min(1, wt/(wd(k) * options.dyxcutoff));
            wt = (cos(wt/max(wt(:))*pi)+1).^2/4;
            tmpdyx = sum(([yy(:)-yyt, xx(:)-xxt]).*wt, 1)./sum(wt);
            dyxt(t,:) = tmpdyx/sum(tmpdyx(:).^2).^0.5;
        end
        yx{k} = yxt(tokeep,:);
        dyx{k} = dyxt(tokeep,:);
        id{k} = k*ones(sum(tokeep(:)), 1, 'single');
    end
    yx = vertcat(yx{:});
    dyx = vertcat(dyx{:});
    id = vertcat(id{:});
    endpt = table(yx, dyx, id);
end