function cropbox = draw_cropbox(IMG)
    hfig = figure('units','normalized','outerposition',[0 0 1 1]);
    imagesc(IMG)
    colormap(gray)
    h = imrect;
    pos = wait(h);
    try close(hfig); end %#ok
    if isempty(pos)
        cropbox = [];
    else
        cropbox = round([pos(2), pos(2)+pos(4), pos(1), pos(1)+pos(3)]);
    end
end

