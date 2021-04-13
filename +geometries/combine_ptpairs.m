function ptpairs = combine_ptpairs(ptpairs1, ptpairs2)
    ptpairs = [ptpairs1; ptpairs2];
    ptpairs = sortrows(ptpairs, 'conf', 'descend');
    Tpt = ptpairs(:, {'yx1', 'yx2'});
    [~, ia, ~] = unique(Tpt, 'first');
    ptpairs = ptpairs(ia, :);
end
