function M = gen_eqtriang_mesh_mask(mask, sideL, extend)
    if nargin < 3
        extend = 0;
    end
    theta = 0;
    imgsz = size(mask);
    imght = imgsz(1);
    imgwd = imgsz(2); 
    vec1 = [sind(theta), cosd(theta)] * sideL;
    vec2 = [sind(theta + 60), cosd(theta + 60)] * sideL;
    vec1n = vec1 - vec2 * dot(vec1, vec2) / dot(vec2, vec2);
    vec2n = vec2 - vec1 * dot(vec1, vec2) / dot(vec1, vec1);
    sideLn2 = dot(vec1n, vec1n);
    gridht0 = imght;
    gridwd0 = imgwd;
    if extend > 0
        gridht = gridht0 + 2 * sideL * (extend + 0.25);
        gridwd = gridwd0 + 2 * sideL * (extend + 0.25);
        if range(mask(:)) > 0
            mask = imdilate(mask, strel('disk', round(sideL * extend)));
        end
    else
        gridht = gridht0;
        gridwd = gridwd0;
    end
    corners = [0, 0; gridht, 0; 0, gridwd; gridht, gridwd];
    cornercoord = corners * [vec1n(:), vec2n(:)] / sideLn2;
    coordrange = ceil(range(cornercoord, 1));
    [vv1, vv2] = meshgrid(1:coordrange(1), 1:coordrange(2));
    V = [vv1(:),vv2(:)];

    pts = V * [vec1; vec2];
    pts = pts - mean(pts, 1) + [gridht0, gridwd0]/2;
    sbs = min(max(1, round(pts)), imgsz);
    ind = sub2ind(imgsz, sbs(:,1), sbs(:,2));
    idx0 = (pts(:, 1) <= (imght + sideL * extend)) &  (pts(:, 1) >= -sideL * extend) & ...
        (pts(:, 2) <= imgwd + sideL * extend) &  (pts(:, 2) >= -sideL * extend);
    idx1 = mask(ind) > 0;
    idx = idx0 & idx1;
    pts0 = pts(idx, :);
    V0 = V(idx,:);
    pts_outside = pts(~idx, :);
    V_outside = V(~idx, :);
    dis = pdist2(pts0, pts_outside, 'euclidean', 'Smallest', 1);
    idxt = dis < sideL * (extend + 0.1);
    yx = [pts0; pts_outside(idxt,:)];
    V = [V0; V_outside(idxt, :)];

    locb0 = (1:size(V,1))';
    [~, locb1] = ismember(V + [0, 1], V, 'rows');
    [~, locb2] = ismember(V + [1, 0], V, 'rows');
    T1 = [locb0, locb1, locb2];
    idxt = all(T1 > 0, 2);
    T1 = T1(idxt, :);
    [~, locb1] = ismember(V + [0, -1], V, 'rows');
    [~, locb2] = ismember(V + [-1, 0], V, 'rows');
    T2 = [locb0, locb1, locb2];
    idxt = all(T2 > 0, 2);
    T2 = T2(idxt, :);
    M = struct;
    M.sideL = sideL;
    M.cList = [T1; T2];
    M.yx = yx;
    M.TR = triangulation([T1;T2], yx(:,2), yx(:,1));
    % M.pt_dis = squareform(pdist(yx)) / sideL;
    M.edges = edges(M.TR);
end