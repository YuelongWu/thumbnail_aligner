function tile = block_tform_render_tile(IMG, MC, row, col, scl, method, export_empty, fillval)
    if nargin < 8
        fillval = 0;
    end
    if nargin < 7
        export_empty = true;
    end
    if nargin < 6
        method = 'linear';
    end
    if nargin < 5
        scl = 1;
    end

    bsz = round(MC.blocksz * scl);
    idxt = find(all(MC.block_id == [row, col], 2), 1);
    if isempty(idxt) || any(isnan(MC.offsets(idxt, :)))
        if export_empty
            tile = zeros(bsz, bsz, 'uint8');
        else
            tile = [];
        end
    else
        [xx, yy] = meshgrid(1:bsz, 1:bsz);
        xy0 = cat(3, xx, yy);
        offset = MC.offsets(idxt,:) * scl;
        deftile = double(MC.def_tiles{idxt});
        if scl ~= 1
            deftile = scl * imresize(deftile, scl);
        end
        deftile = deftile + reshape(offset,1,1,2);
        deftile_flat = reshape(deftile, size(deftile,1) * size(deftile, 2), 2);
        minxy = min(deftile_flat, [], 1);
        maxxy = max(deftile_flat, [], 1);
        minxy = max(floor(minxy) - 1, 1);
        maxxy = min(ceil(maxxy) + 1, [size(IMG,2), size(IMG,1)]);
        IMGt = IMG(minxy(2):maxxy(2), minxy(1):maxxy(1));
        if isempty(IMGt)
            tile = [];
            return;
        end
        deftile = deftile - reshape(minxy,1,1,2) + 1;
        D = deftile - xy0;
        IMGt(IMGt == fillval) = IMGt(IMGt == fillval) + sign(127 - fillval);
        tile = imwarp(IMGt, D, method, 'FillValues', fillval);
    end
end