function [ptpairs, modified] = edit_ptpairs(IMG1, IMG2, ptpairs, masks, cropbox, catdim)
    rdg = false;
    modified = false;
    showconf = false;
    if isstruct(IMG1)
        IMG1 = IMG1.img;
    end
    if isstruct(IMG2)
        IMG2 = IMG2.img;
    end
    if nargin < 5 || isempty(cropbox)
        cropbox = [1, 1, inf, inf];
    end
    IMG1 = local_crop_img(IMG1, cropbox);
    IMG2 = local_crop_img(IMG2, cropbox);
    crpshft = cropbox(1:2) - 1;
    if ~isempty(ptpairs)
        ptpairs = local_shift_ptpairs(ptpairs, -crpshft);
    end
    screen_ratio = 1.6;
    imght1 = size(IMG1, 1);
    imgwd1 = size(IMG1, 2);
    imght2 = size(IMG2, 1);
    imgwd2 = size(IMG2, 2);
    if nargin < 6 || isempty(catdim)
        ratio1 = max(imgwd1, imgwd2)/(imght1 + imght2);
        ratio2 = (imgwd1 + imgwd2)/max(imght1, imght2);
        [~, catdim] = min([abs(log(ratio1/screen_ratio)), abs(log(ratio2/screen_ratio))]);
    end


    if nargin < 4 || isempty(masks)
        masks = [];
    else
        if ischar(masks{1})
            masks{1} = imread(masks{1});
        end
        if ischar(masks{2})
            masks{2} = imread(masks{2});
        end
        masks{1} = local_crop_img(masks{1}, cropbox);
        masks{2} = local_crop_img(masks{2}, cropbox);
        IMG1(masks{1} == 0) = 0;
        IMG2(masks{2} == 0) = 0;
        masks{1} = mask_gen.expand_mask(masks{1});
        masks{2} = mask_gen.expand_mask(masks{2});
    end
    if catdim == 1
        imgwd = max(imgwd1, imgwd2);
        IMG = [padarray(IMG1, [0, imgwd-imgwd1], 'post'); ...
            padarray(IMG2, [0, imgwd-imgwd2], 'post')];
        dyx2 = [imght1, 0];
    else
        imght = max(imght1, imght2);
        IMG = [padarray(IMG1, [imght-imght1,0], 'post'), ...
            padarray(IMG2, [imght-imght2,0], 'post')];
        dyx2 = [0, imgwd1];
    end
    hfig = figure(923);
    shline = false;
    finished = false;
    blk = 0;
    firsttime = true;
    while ~finished
        switch blk
        case 2
            hfig = visualize_ptpairs(hfig, IMG, ptpairs, dyx2, shline);
        case 1
            hfig = figure(923); local_visualize_matching_b(hfig, IMG1, ptpairs, false);
        case 0
            hfig = figure(923); local_visualize_matching(hfig, IMG1, IMG2, ptpairs, false);
        end
        if showconf
            text(double(ptpairs.yx1(:,2)), double(ptpairs.yx1(:,1)), cellstr(num2str(ptpairs.conf,'%.2f')), 'color','w');
        end
        % return
        if firsttime
            disp('>>');
            firsttime = false;
        end
        try
            w = waitforbuttonpress;
        catch ME
            if strcmpi(ME.identifier, 'MATLAB:UI:CancelWaitForKeyOrButtonPress')
                modified = -1;
                return
            end
        end
        try
            if w
                switch lower(get(hfig,'CurrentCharacter'))
                case 'q'    % exit
                    ptpairs = local_shift_ptpairs(ptpairs, crpshft);
                    if ~isempty(ptpairs)
                        ptpairs.conf = ones(size(ptpairs, 1),1);
                    end
                    disp('exit manual mode...')
                    finished = true;
                case 'u'    % update visualization
                    blk = blk + 1;
                    if blk >= 3
                        blk = blk - 3;
                    end
                    % hfig = visualize_ptpairs(hfig, IMG, ptpairs, dyx2, shline);
                case 'a'    % addpoint
                    ptpairs = add_points(hfig, ptpairs, dyx2, masks);
                    modified = true;
                case 'r'    % remove points
                    ptpairs = remove_points(hfig, ptpairs, dyx2);
                    modified = true;
                case 'c'    % clear points
                    flag = input('sure to clear?(Y/N): ', 's');
                    if strcmpi(flag, 'y')
                        ptpairs = ptpairs([], :);
                    end
                    modified = true;
                case 'l'
                    shline = ~shline;
                case 'f'
                    hfig2 = figure('units','normalized','outerposition',[0.05 0.05 0.9 0.9]);
                    local_visualize_matching(hfig2, IMG1, IMG2, ptpairs, rdg);
                    uiwait(hfig2);
                case 'b'
                    hfig2 = figure('units','normalized','outerposition',[0.05 0.05 0.9 0.9]);
                    local_visualize_matching_b(hfig2, IMG1, ptpairs, rdg);
                    uiwait(hfig2)
                case 'y'
                    d = input('# of pairs to delete: ');
                    if ~isempty(d)
                        A = geometries.fit_affine(ptpairs.yx1, ptpairs.yx2);
                        yx2t = ptpairs.yx2 * A(1:2,1:2) + A(3,1:2);
                        dis = sum((ptpairs.yx1 - yx2t).^2, 2);
                        [~, idx] = maxk(dis, d);
                        ptpairs(idx, :) = [];
                    end
                    modified = true;
                case 'o'
                    ptpairs = geometries.filter_iter_outlier(ptpairs);
                    modified = true;
                case 'p'
                    ptpairs = geometries.filter_fit_surface(ptpairs);
                    modified = true;
                case 'm'
                    flag = input('sure to delete file?', 's');
                    if strcmpi(flag, 'y')
                        modified = -2;
                        return
                    end
                otherwise
                end
            end
        catch ME
            disp(ME.message);
        end
    end
end

function ptpairs = add_points(hfig, ptpairs, dyx2, masks)
    if ~isvalid(hfig) || isempty(hfig.Children)
        return
    end
    [xi, yi] = getpts(hfig);
    yxi = [yi(:), xi(:)];
    idx2 = all((yxi - dyx2) > 0, 2);
    yx2 = yxi(idx2, :) - dyx2;
    yx1 = yxi(~idx2, :);
    if size(yx1, 1) ~= size(yx2, 1)
        return
    end
    N = size(yx1, 1);
    if isempty(masks)
        region_id = ones(N, 2, 'single');
    else
        mask = masks{1};
        rid1 = get_region_id(yx1, mask);
        mask = masks{2};
        rid2 = get_region_id(yx2, mask);
        region_id = [rid1(:), rid2(:)];
    end
    conf = 1 * ones(N, 1, 'single');
    rotation = nan(N, 1, 'single');
    if isempty(ptpairs)
        ptpairs = table(yx1, yx2, region_id, conf, rotation);
    else
        ptpairs = [table(yx1, yx2, region_id, conf, rotation); ptpairs];
    end
end

function ptpairs = remove_points(hfig, ptpairs, dyx2)
    if ~isvalid(hfig) || isempty(hfig.Children) || isempty(ptpairs)
        return
    end
    ax = hfig.Children(1);
    yx1 = ptpairs.yx1;
    yx2 = ptpairs.yx2 + dyx2;
    h = imrect(ax); position = wait(h); delete(h);
    crn1 = [position(2), position(1)];
    crn2 = [position(4), position(3)] + crn1;
    idx1 = all((yx1 >= crn1) & (yx1 < crn2),2);
    idx2 = all((yx2 >= crn1) & (yx2 < crn2),2);
    idx = idx1 | idx2;
    ptpairs = ptpairs(~idx,:);
end

function hfig = visualize_ptpairs(hfig, IMG, ptpairs, dyx2, shline)
    if isvalid(hfig)
        clf(hfig);
        ax = axes(hfig);
    else
        hfig = figure;
        ax = axes(hfig);
    end
    imagesc(ax, IMG);
    set(ax, 'Position',[0,0,1,1]);
    colormap(hfig, gray);
    if isempty(ptpairs)
        return
    end
    yx1 = ptpairs.yx1;
    yx2 = ptpairs.yx2 + dyx2;
    [shft, dim] = max(dyx2);
    idx1 = yx1(:,dim) >= 0 & yx1(:,dim) <= shft;
    idx2 = yx2(:,dim) >= shft & yx1(:,dim) <= size(IMG, dim);
    idx = idx1 & idx2;
    xx = [yx1(idx,2), yx2(idx,2)]';
    yy = [yx1(idx,1), yx2(idx,1)]';
    hold on;
    if shline
        plot(ax, xx, yy, 'b-');
    end
    plot(ax, yx1(idx,2), yx1(idx,1), 'r.','MarkerSize',10);
    plot(ax, yx2(idx,2), yx2(idx,1), 'g.','MarkerSize',10);
end

function IMG2t = local_visualize_matching(hfig, IMG1, IMG2, ptpairs, rgd)
    if nargin < 5
        rgd = false;
    end
    xy1 = fliplr(ptpairs.yx1);
    xy2 = fliplr(ptpairs.yx2);
    if rgd
        [~, R] = geometries.fit_affine(xy1, xy2);
    else
        [R, ~] = geometries.fit_affine(xy1, xy2);
    end
    xy2t = xy2 * R(1:2,1:2) + R(3,1:2);
    IMG2t = imwarp(IMG2, affine2d(R), 'OutputView', imref2d(size(IMG1)));
    clf;
    ax = axes(hfig);
    % imagesc(ax, (255-histeq(cat(3, IMG1, IMG2t, IMG1))));
    imagesc(ax, 255 - 2*(255-(cat(3, IMG1, IMG2t, IMG1))));
    set(ax, 'Position',[0,0,1,1]);
    hold on;
    axis off
    xx = [xy1(:,1), xy2t(:,1)]';
    yy = [xy1(:,2), xy2t(:,2)]';
    hold on;

    plot(ax, xx, yy, 'y-');
    plot(ax, xy1(:,1), xy1(:,2), 'ro', 'MarkerSize',3);
    plot(ax, xy2t(:,1), xy2t(:,2), 'go','MarkerSize',5);
end

function local_visualize_matching_b(hfig, IMG1, ptpairs, rgd)
    if nargin < 5
        rgd = false;
    end
    xy1 = fliplr(ptpairs.yx1);
    xy2 = fliplr(ptpairs.yx2);
    if rgd
        [~, R] = geometries.fit_affine(xy1, xy2);
    else
        [R, ~] = geometries.fit_affine(xy1, xy2);
    end
    xy2t = xy2 * R(1:2,1:2) + R(3,1:2);

    clf(hfig)
    ax = axes(hfig);
    image(ax, 0*IMG1);
    colormap('gray')
    set(ax, 'Position',[0.05,0.05,0.9,0.9]);
    set(ax, 'Color',[0, 0, 0]);
    hold on;
    xx = [xy1(:,1), xy2t(:,1)]';
    yy = [xy1(:,2), xy2t(:,2)]';
   
    plot(ax, xx, yy, 'y-');
    plot(ax, xy1(:,1), xy1(:,2), 'r.');
    plot(ax, xy2t(:,1), xy2t(:,2), 'g.');

    %
    % dis = sum((xy2t - xy1).^2,2).^0.5;
    % cc = (dis - mean(dis)) / std(dis);
    % text(double(xy1(:,1)), double(xy1(:,2)), cellstr(num2str(cc(:),'%.2f')), 'color','w')

    axis off
end

function rids = get_region_id(yx, mask)
    imght = size(mask, 1);
    imgwd = size(mask, 2);
    yy = min(imght, max(1, round(yx(:, 1))));
    xx = min(imgwd, max(1, round(yx(:, 2))));
    ind = sub2ind(size(mask), yy, xx);
    rids = mask(ind);
end

function crpimg = local_crop_img(img, cropbox)
    cropbox(1) = max(1, cropbox(1));
    cropbox(2) = max(1, cropbox(2));
    cropbox(3) = min(size(img,1), cropbox(3));
    cropbox(4) = min(size(img,2), cropbox(4));
    crpimg = img(cropbox(1):cropbox(3), cropbox(2):cropbox(4));
end

function ptpairs = local_shift_ptpairs(ptpairs, dyx)
    ptpairs.yx1 = ptpairs.yx1 + dyx;
    ptpairs.yx2 = ptpairs.yx2 + dyx;
end