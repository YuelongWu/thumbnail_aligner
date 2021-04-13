function indx = add_lines_enpts(imgsz, start_yx, end_yx, return_matrix, shoulder_flag, flipin)
    if nargin < 6
        flipin = false;
    end
    if nargin < 5
        shoulder_flag = true;
    end
    if nargin < 4
        return_matrix = false;
    end
    if flipin
        start_yx = fliplr(start_yx);
        end_yx = fliplr(end_yx);
    end
    Np = size(start_yx, 1);
    indx = cell(Np, 1);
    imght = imgsz(1);
    imgwd = imgsz(2);
    start_yx = round(start_yx);
    end_yx = round(end_yx);
    maxlen = 0;
    entr1 = nan(Np, 1);
    for k = 1:Np
        x0 = start_yx(k,2);
        y0 = start_yx(k,1);
        x1 = end_yx(k,2);
        y1 = end_yx(k,1);
        dx = sign(x1 - x0);
        dy = sign(y1 - y0);
        Nl = max(abs(y1 - y0), abs(x1 - x0)) + 1;
        xx = round(linspace(x0, x1, Nl))';
        yy = round(linspace(y0, y1, Nl))';
        if shoulder_flag
            jump_idx = abs(diff(xx)) > 0 & abs(diff(yy)) > 0;
            if any(jump_idx)
                xi = xx(jump_idx);
                yi = yy(jump_idx);
                mask = mod(1:numel(xi), 2); % alternating shoulder points
                xi = xi(:) + dx * mask(:);
                yi = yi(:) + dy * (1 - mask(:));
                xx = [xx(:); xi(:)];
                yy = [yy(:); yi(:)];
            end
        end
        idx = (xx > 0) & (yy > 0) & (xx <= imgwd) & (yy <= imght);
        xx = xx(idx);
        yy = yy(idx);
        inds = sub2ind(imgsz, yy, xx);
        if isempty(inds)
            continue
        end
        entr1(k) = inds(1);
        maxlen = max(numel(inds), maxlen);
        indx{k} = inds;
    end
    if return_matrix
        indx0 = indx;
        indx = repmat(entr1, 1, maxlen);
        for k = 1:Np
            inds = indx0{k};
            if ~isempty(inds)
                indx(k, 1:numel(inds)) = inds;
            end
        end
    end
end
