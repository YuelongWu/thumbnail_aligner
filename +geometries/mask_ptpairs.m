function ptpairs = mask_ptpairs(ptpairs, mask1, mask2)
    yx1 = ptpairs.yx1;
    yx2 = ptpairs.yx2;
    idx = true(size(yx1,1),1);
    if ~isempty(mask1)
        imgsz = [size(mask1,1), size(mask1,2)];
        yx1 = round(yx1);
        idx = idx & all(yx1>0 & yx1 < imgsz,2);
        yx1 = min(max(yx1,1),imgsz);
        ind = sub2ind(imgsz, yx1(:,1), yx1(:,2));
        idx = idx & (mask1(ind) > 0);
    end
    if ~isempty(mask2)
        imgsz = [size(mask2,1), size(mask2,2)];
        yx2 = round(yx2);
        idx = idx & all(yx2>0 & yx2 < imgsz,2);
        yx2 = min(max(yx2,1),imgsz);
        ind = sub2ind(imgsz, yx2(:,1), yx2(:,2));
        idx = idx & (mask2(ind) > 0);
    end
    ptpairs = ptpairs(idx,:);
end