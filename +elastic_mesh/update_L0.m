function M = update_L0(M)
    xx = M.TR.Points(:, 1);
    yy = M.TR.Points(:, 2);
    edg_idx = M.edges;
    dx = diff(xx(edg_idx),1,2);
    dy = diff(yy(edg_idx),1,2);
    L = (dx.^2 + dy.^2).^0.5;
    M.L0 = L;
end