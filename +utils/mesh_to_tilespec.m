function sout = mesh_to_tilespec(Ms, imgpaths, imgsz, layer, subsample_ratio)
    if nargin < 5 || isempty(subsample_ratio)
        subsample_ratio = 1;
    end
    if nargin < 4 || isempty(layer)
        layer = 1;
    end
    if nargin < 3 || isempty(imgsz)
        IMG0 = imread(imgpaths{1});
        imgsz = size(IMG0);
    end

    s0 = template_str(false);
    s0 = strrep(s0, '$HEIGHT$', num2str(imgsz(1)));
    s0 = strrep(s0, '$WIDTH$', num2str(imgsz(2)));
    s0 = strrep(s0, '$LAYER$', num2str(layer));

    S = cell(numel(Ms), 1);
    offset = [nan, nan];
    for k = 1 : numel(Ms)
        M = Ms{k};
        xy = min(M.TR.Points, [], 1);
        offset = min(offset, xy);
    end
    for k = 1 : numel(Ms)
        if k < numel(Ms)
            s = s0;
        else
            s = template_str(true);
            s = strrep(s, '$HEIGHT$', num2str(imgsz(1)));
            s = strrep(s, '$WIDTH$', num2str(imgsz(2)));
            s = strrep(s, '$LAYER$', num2str(layer));
        end
        s = strrep(s, '$FILE$', imgpaths{k});
        M = Ms{k};
        xy0 = M.TR0.Points;
        xy1 = M.TR.Points - offset;
        xy_min = round(min(xy1, [], 1), 1);
        xy_max = round(max(xy1, [], 1), 1);
        s = strrep(s, '$X0$', num2str(xy_min(1)));
        s = strrep(s, '$X1$', num2str(xy_max(1)));
        s = strrep(s, '$Y0$', num2str(xy_min(2)));
        s = strrep(s, '$Y1$', num2str(xy_max(2)));
        s = strrep(s, '$TILE_INDEX$', num2str(k));
        DXY = mean(xy1, 1) - mean(xy0, 1);
        XY = [xy0, xy1 - DXY, ones(size(xy0,1),1)];
        if subsample_ratio < 1
            ctr_xy = fliplr(imgsz) / 2;
            xyr = 1 - 2 * max(abs(xy0(:,1) - ctr_xy(:,1))/imgsz(2), ...
                              abs(xy0(:,2) - ctr_xy(:,2))/imgsz(1));
            xyr = max(0, xyr - 0.15);
            if max(xyr(:)) > 0
                xyr = xyr / max(xyr(:)) * (1-subsample_ratio);
                rndr = rand(size(xyr));
                idxt = rndr > xyr;
                XY = XY(idxt,:);
            end
        end
        XY = XY.';
        datastrC = erase(cellstr(num2str(XY(:),'%.1f')),' ');
        datastr = strjoin(datastrC, ' ');
        s = strrep(s, '$PTS$', datastr);
        s = strrep(s, '$DX$', num2str(DXY(1), '%.1f'));
        s = strrep(s, '$DY$', num2str(DXY(2), '%.1f'));
        S{k} = s;
    end
    sout = ['[', newline, horzcat(S{:}), ']', newline];
end

function s = template_str(last_line)
    tb = '    ';
    nl = newline;
    s = [tb, '{', nl, ...
        tb, tb, '"bbox": [', nl, ...
        tb, tb, tb, '$X0$,', nl, ...
        tb, tb, tb, '$X1$,', nl, ...
        tb, tb, tb, '$Y0$,', nl, ...
        tb, tb, tb, '$Y1$', nl, ...
        tb, tb, '],', nl, ...
        tb, tb, '"height": $HEIGHT$,', nl, ...
        tb, tb, '"layer": $LAYER$,', nl, ...
        tb, tb, '"maxIntensity": 255.0,', nl, ...
        tb, tb, '"minIntensity": 0.0,', nl, ...
        tb, tb, '"mipmapLevels": {', nl, ...
        tb, tb, tb, '"0": {', nl, ...
        tb, tb, tb, tb, ' "imageUrl": "file://$FILE$"', nl, ...
        tb, tb, tb, '}', nl, ...
        tb, tb, '},', nl, ...
        tb, tb, '"tile_index": $TILE_INDEX$,', nl, ...
        tb, tb, '"transforms": [', nl, ...
        tb, tb, tb, '{', nl, ...
        tb, tb, tb, tb, '"className": "mpicbg.trakem2.transform.RigidModel2D",', newline, ...
        tb, tb, tb, tb, '"dataString": "0 $DX$ $DY$"', newline, ...
        tb, tb, tb, '},', newline, ...
        tb, tb, tb, '{', nl, ...
        tb, tb, tb, tb, '"className": "mpicbg.trakem2.transform.PointsTransformModel",', newline, ...
        tb, tb, tb, tb, '"dataString": "$PTS$"', newline, ...
        tb, tb, tb, '}', newline, ...
        tb, tb, '],', newline, ...
        tb, tb, '"width": $WIDTH$', newline];
        if last_line
            s = [s, tb, '}', newline];
        else
            s = [s, tb, '},', newline];
        end
end