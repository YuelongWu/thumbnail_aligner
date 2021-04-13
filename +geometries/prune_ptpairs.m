function ptpairs = prune_ptpairs(ptpairs)
    ptpairs = sortrows(ptpairs, 'conf', 'descend');
    Tpt = ptpairs(:, 'yx1');
    [~, ia, ~] = unique(Tpt, 'first');
    ptpairs = ptpairs(ia, :);
    Tpt = ptpairs(:, 'yx2');
    [~, ia, ~] = unique(Tpt, 'first');
    ptpairs = ptpairs(ia, :);
end