function D = geodesic_dis_diff(mask, seed_mask, scl)
    if nargin < 3
        scl = 1;
    end
    masksz = size(mask);
    if scl ~= 1
        mask = imresize(single(mask), scl, 'bilinear') > 0;
        seed_mask = imresize(single(seed_mask), scl, 'bilinear') > 0;
    end
    Dis1 = bwdistgeodesic(mask, seed_mask, 'quasi-euclidean');
    Dis2 = bwdist(seed_mask, 'quasi-euclidean');
    Dis_diff = Dis1 - Dis2;
    D = imresize(Dis_diff, masksz, 'nearest') / scl;
end
