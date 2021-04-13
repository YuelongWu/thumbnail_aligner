function ptpairs = flip_ptpairs(ptpairs)
    if isempty(ptpairs)
        return
    end
    yx2 = ptpairs.yx1;
    ptpairs.yx1 = ptpairs.yx2;
    ptpairs.yx2 = yx2;
    ptpairs.region_id = fliplr(ptpairs.region_id);
end
