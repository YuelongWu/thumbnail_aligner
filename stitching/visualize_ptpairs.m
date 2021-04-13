function visualize_ptpairs(PTPAIRS, Ms,imgpaths)
    As = zeros(3,3, numel(Ms));
    for k = 1:numel(Ms)
        As(:,:,k) = geometries.fit_affine(Ms{k}.TR.Points, Ms{k}.TR0.Points);
    end
    % clf;
    hold on
    ds = 1;
    for k = 1:numel(PTPAIRS)
        mnames = PTPAIRS(k).imgnames;
        xy1 = fliplr(PTPAIRS(k).ptpairs.yx1);
        xy2 = fliplr(PTPAIRS(k).ptpairs.yx2);
        idx1 = find(contains(imgpaths, mnames{1}),1);
        idx2 = find(contains(imgpaths, mnames{2}),1);
        xy1t = xy1 * As(1:2,1:2,idx1) + As(3,1:2, idx1);
        xy2t = xy2 * As(1:2,1:2,idx2) + As(3,1:2, idx2);
        cc = hsv2rgb([rand(1), 0.25 + 0.5*rand(1,2)]);
        plot([xy1t(1:ds:end,1), xy2t(1:ds:end,1)]',[xy1t(1:ds:end,2), xy2t(1:ds:end,2)]', 'Color', cc)
        plot(xy1t(1:ds:end,1),xy1t(1:ds:end,2),'.', 'Color', cc)
    end
end