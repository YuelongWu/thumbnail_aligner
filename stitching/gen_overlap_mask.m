function [mask1, mask2] = gen_overlap_mask(yx1, yx2, imgsz1, imgsz2, options)
    if nargin < 5
        options = [];
    end
    if nargin < 4 || isempty(imgsz2)
        imgsz2 = imgsz1;
    end
    options = utils.set_default(options, 'margin', 0);
    options = utils.set_default(options, 'whisker', 75);
    yx2 = round(yx2 - yx1 + 1);
    yx1 = ones(size(yx1));

    ht1 = imgsz1(1);
    wd1 = imgsz1(2);
    ht2 = imgsz2(1);
    wd2 = imgsz2(2);

    mask1 = false(ht1, wd1);
    mask2 = false(ht2, wd2);

    mrgn = options.margin;
    if any([ht1, wd1, ht2, wd2] <= 2*mrgn)
        return
    end
    % bbox: [ymin, ymax, xmin, xmax]
    mrgn_step = mrgn * [1, -1, 1, -1];
    bbox1 = [yx1(1), yx1(1) + ht1 - 1, yx1(2), yx1(2) + wd1 - 1] + mrgn_step;
    bbox2 = [yx2(1), yx2(1) + ht2 - 1, yx2(2), yx2(2) + wd2 - 1] + mrgn_step;
    [bbox, ovlp] = rect_intersect(bbox1, bbox2);
    if ovlp == 0
        return
    end
    whisker = options.whisker;
    whisker_step = whisker * [-1, 1, -1, 1];
    bbox_t = bbox + whisker_step;
    idxt1 = rect_intersect(bbox1, bbox_t);
    idxt2 = rect_intersect(bbox2, bbox_t);
    idxt1 = idxt1(:) - repelem(yx1(:), 2) + 1;
    idxt2 = idxt2(:) - repelem(yx2(:), 2) + 1;
    mask1(idxt1(1):idxt1(2), idxt1(3):idxt1(4)) = 1;
    mask2(idxt2(1):idxt2(2), idxt2(3):idxt2(4)) = 1;
end
