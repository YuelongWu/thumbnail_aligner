function IDs = build_2d_tree(yx, ids, leaf_size, range_yx)
    yx_std = nanstd(yx, 1, 1);
    if abs(numel(ids)/2 - leaf_size) >= abs(numel(ids) - leaf_size) || max(yx_std) <= range_yx
        IDs = {ids(:)};
        return
    end
    [~, longer_dim] = max(yx_std);
    idx1 = yx(:, longer_dim) >= median(yx(:, longer_dim), 'omitnan');
    idx2 = yx(:, longer_dim) < median(yx(:, longer_dim), 'omitnan');
    if all(idx1) || all(idx2)
        IDs = {ids(:)};
        return
    end
    IDs1 = geometries.build_2d_tree(yx(idx1,:), ids(idx1), leaf_size, range_yx);
    IDs2 = geometries.build_2d_tree(yx(idx2,:), ids(idx2), leaf_size, range_yx);
    IDs = [IDs1; IDs2];
end
