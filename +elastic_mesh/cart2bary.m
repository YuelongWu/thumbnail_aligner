function [B, ID, oob] = cart2bary(TR, yx)
    xy = double(fliplr(yx));
    [ID, B] = pointLocation(TR,xy);
    oob = isnan(ID);
    if any(oob)
        c_xy = incenter(TR);
        xy_t = xy(oob, :);
        dis = pdist2(xy_t, c_xy);
        [~, min_idx] = min(dis, [], 2);
        ID(oob) = min_idx;
        B(oob, :) = cartesianToBarycentric(TR,min_idx,xy_t);
    end
end