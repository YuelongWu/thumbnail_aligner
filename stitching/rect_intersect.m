function [bbox, itsect] = rect_intersect(bbox1, bbox2)
    % bbox: ymin, ymax, xmin, xmax
    idx1 = max(bbox1(:,1), bbox2(:,1));
    idx2 = min(bbox1(:,2), bbox2(:,2));
    idx3 = max(bbox1(:,3), bbox2(:,3));
    idx4 = min(bbox1(:,4), bbox2(:,4));
    bbox = [idx1, idx2, idx3, idx4];
    itsect = max(0, idx2 - idx1) .* max(0, idx4 - idx3);
end