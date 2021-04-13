function [PTPAIRS, modified] = edit_ptpairs_stitching(PTPAIRS, imgpath, imgsz0, options)
    if nargin < 4 || isempty(options)
        options = struct;
    end
    if nargin < 3
        imgsz0 = [];
    end
    options = utils.set_default(options, 'err_thresh', 0);
    options = utils.set_default(options, 'scale', 0.25);
    modified = 0;
    kept = true(numel(PTPAIRS), 1);
    for k = 1:numel(PTPAIRS)
        ptpr = PTPAIRS(k).ptpairs;
        if isempty(ptpr)
            kept(k) = false;
            continue
        end
        yx1 = ptpr.yx1;
        yx2 = ptpr.yx2;
        [~, R] = geometries.fit_affine(yx1, yx2);
        yx2t = yx2 * R(1:2,1:2) + R(3,1:2);
        dis = sum((yx1 - yx2t).^2, 2).^0.5;
        if max(dis(:)) <= options.err_thresh
            continue
        end
        IMG1 = [];
        IMG2 = [];
        if ~isempty(imgpath)
            if exist([imgpath, filesep, PTPAIRS(k).imgnames{1}], 'file')
                IMG1 = imread([imgpath, filesep, PTPAIRS(k).imgnames{1}]);
                if options.scale ~= 1
                    IMG1 = imresize(IMG1, options.scale);
                end
            end
            if exist([imgpath, filesep, PTPAIRS(k).imgnames{2}], 'file')
                IMG2 = imread([imgpath, filesep, PTPAIRS(k).imgnames{2}]);
                if options.scale ~= 1
                    IMG2 = imresize(IMG2, options.scale);
                end
            end
        end
        if isempty(imgsz0)
            imgsz = ceil(max(max(yx1, [], 1), max(yx2, [], 1)));
            imgsz = max(ceil(1.05 * imgsz), imgsz + 50);
        else
            imgsz = imgsz0;
        end
        if options.scale ~= 1
            imgsz = round(options.scale * imgsz);
        end
        if isempty(IMG1)
            continue
            % IMG1 = 50 * ones(imgsz, 'uint8');
        end
        if isempty(IMG2)
            continue
            % IMG2 = 50 * ones(imgsz, 'uint8');
        end
        disp([PTPAIRS(k).imgnames{1}, ' <-> ',PTPAIRS(k).imgnames{2}])
        if options.scale ~= 1
            ptpr.yx1 = ptpr.yx1 * options.scale;
            ptpr.yx2 = ptpr.yx2 * options.scale;
        end
        [ptpr, mdfed] = edit_ptpairs(IMG1, IMG2, ptpr);
        if options.scale ~= 1
            ptpr.yx1 = ptpr.yx1 / options.scale;
            ptpr.yx2 = ptpr.yx2 / options.scale;
        end
        if mdfed == -1
            modified = -1;
            return
        elseif mdfed
            if isempty(ptpr)
                kept(k) = false;
            else
                PTPAIRS(k).ptpairs = ptpr;
            end
            modified = 1;
        end
    end
    PTPAIRS = PTPAIRS(kept);
end
