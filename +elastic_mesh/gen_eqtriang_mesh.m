function M = gen_eqtriang_mesh(imgsz, sideL, full_cover, return_mesh)
    if nargin < 4
        return_mesh = true;
    end
    if nargin < 3
        full_cover = false;
    end
    theta = 0;
    imght = imgsz(1);
    imgwd = imgsz(2); 
    vec1 = [sind(theta), cosd(theta)] * sideL;
    vec2 = [sind(theta + 60), cosd(theta + 60)] * sideL;
    
    vec1n = vec1 - vec2 * dot(vec1, vec2) / dot(vec2, vec2);
    vec2n = vec2 - vec1 * dot(vec1, vec2) / dot(vec1, vec1);
    sideLn2 = dot(vec1n, vec1n);
    gridht0 = imght;
    gridwd0 = imgwd;
    if full_cover
        gridht = gridht0 + 2 * sideL;
        gridwd = gridwd0 + 2 * sideL;
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
    pts = pts - mean(pts, 1) + [gridht, gridwd]/2;
    idx = (pts(:, 1) <= gridht) &  (pts(:, 1) >= 0) & (pts(:, 2) <= gridwd) &  (pts(:, 2) >= 0);
    pts = pts(idx, :);
    yx = pts - [gridht, gridwd]/2 + [imght, imgwd]/2;
    if ~return_mesh
        M = yx;
        return
    end
    V = V(idx, :);
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
    M.pt_dis = squareform(pdist(yx)) / sideL;
    M.edges = edges(M.TR);
end
