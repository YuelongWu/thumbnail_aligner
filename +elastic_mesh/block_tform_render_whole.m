function [IMGt, Dt] = block_tform_render_whole(IMG, MC, scl, method, fillval)
    if nargin < 5
        fillval = 0;
    end
    if nargin < 4 || isempty(method)
        method = 'linear';
    end
    if nargin < 3
        scl = 1;
    end
    outsz = round(MC.outsz * scl);
    bsz = round(MC.blocksz * scl);
    IMGt = fillval * ones(outsz, class(IMG));
    if nargout > 1
        Dt = nan(outsz(1), outsz(2), 2);
    end
    [xx, yy] = meshgrid(1:bsz, 1:bsz);
    xy0 = cat(3, xx, yy);
    for k = 1:MC.Nblk
        offset = MC.offsets(k,:) * scl;
        if any(isnan(offset))
            continue
        end
        deftile = double(MC.def_tiles{k});
        bx = MC.bbox(k,:);
        if scl ~= 1
            deftile = scl * imresize(deftile, scl);
            bx([1,3]) = round((bx([1,3]) - 1) * scl) + 1;
            bx([2,4]) = round(bx([2,4]) * scl);
        end
        deftile = deftile + reshape(offset,1,1,2);
        D = deftile - xy0;
        tile = imwarp(IMG, D, method, 'FillValues', fillval);
        IMGt(bx(1):bx(2), bx(3):bx(4)) = tile;
        if nargout > 1
            Dt(bx(1):bx(2), bx(3):bx(4),:) = deftile;
        end
    end
    if nargout > 1
        [xx, yy] = meshgrid(1:outsz(2), 1:outsz(1));
        xy0 = cat(3, xx, yy);
        Dt = Dt - xy0;
        Dx = Dt(:,:,1);
        nanmask = isnan(Dx);
        if any(nanmask(:))
            Dy = Dt(:,:,2);
            [~,idxt] = bwdist(~nanmask);
            Dt(:,:,1) = Dx(idxt);
            Dt(:,:,2) = Dy(idxt);
        end
    end
end
