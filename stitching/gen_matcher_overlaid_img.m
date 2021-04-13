function [IMGc, mxdis] = gen_matcher_overlaid_img(IMGt, PTPAIRS, Ms, imglist, scl)
    if numel(scl) == 1
        scl = [scl, scl];
    end
    if isstruct(imglist)
        imgnames = {imglist.name};
        imgnames = imgnames(:);
    else
        imgnames = cell(numel(imglist),1);
        for t = 1:numel(imglist)
            [~, nmt, ext] = fileparts(imglist{t});
            imgnames{t} = [nmt, ext];
        end
    end
    YXs = cell(numel(PTPAIRS), 1);
    mxdis = 0;
    for k = 1:numel(PTPAIRS)
        ptpr = PTPAIRS(k).ptpairs;
        nms = PTPAIRS(k).imgnames;
        idxt = find(strcmpi(imgnames, nms{1}), 1);
        M = Ms{idxt};
        yx1 = ptpr.yx1;
        myx0 = fliplr(M.TR0.Points);
        myx1 = fliplr(M.TR.Points);
        myx1 = myx1(1:size(myx0,1),:);
        A = geometries.fit_affine(myx1, myx0);
        yx1t = yx1 * A(1:2,1:2) + A(3,1:2);
        YXs{k} = yx1t;
        idxt = find(strcmpi(imgnames, nms{2}), 1);
        M = Ms{idxt};
        yx2 = ptpr.yx2;
        myx0 = fliplr(M.TR0.Points);
        myx1 = fliplr(M.TR.Points);
        myx1 = myx1(1:size(myx0,1),:);
        A = geometries.fit_affine(myx1, myx0);
        yx2t = yx2 * A(1:2,1:2) + A(3,1:2);
        dis0 = max(sum((yx1t - yx2t) .^ 2, 2).^0.5);
        if dis0 > mxdis
            mxdis = dis0;
        end
    end
    yxs = vertcat(YXs{:});

    IMGc = 0.75 * imresize(IMGt, scl(1));
    yxs = min(max(1, round(yxs * scl(2))), size(IMGc));
    ind = sub2ind(size(IMGc), yxs(:,1), yxs(:,2));
    IMGct = IMGc;
    IMGct(ind) = 255;
    IMGc = cat(3, IMGct, IMGc, IMGc);
end