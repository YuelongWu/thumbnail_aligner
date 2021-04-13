function xy = bary2cart(TR, ID, B, offst)
    if nargin < 4
        offst = 0;
    end
    cList = TR.ConnectivityList(ID, :) + offst;
    xyt = reshape(TR.Points(cList,:), [size(cList,1),size(cList,2),2]);
    xyn = sum(xyt .* B, 2);
    xy = permute(xyn, [1,3,2]);
end