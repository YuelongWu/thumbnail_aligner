function [Ms, outsz] = reset_mesh_origin(Ms, offset0)
    if nargin < 2
        offset0 = 0;
    end
    xy_max = [nan, nan];
    xy_min = [nan, nan];
    for t = 1:numel(Ms)
        M = Ms{t};
        xy_min = min(xy_min, min(M.TR.Points, [], 1));
        xy_max = max(xy_max, max(M.TR.Points, [], 1));
    end
    offset = offset0 - xy_min;
    outsz = ceil(fliplr(xy_max + offset + offset0));
    for t = 1:numel(Ms)
        Ms{t}.TR.Points = Ms{t}.TR.Points + offset;
    end
end
