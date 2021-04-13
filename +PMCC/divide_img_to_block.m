function [blocks, cntr_yx] = divide_img_to_block(IMG, blocksz, fillval)
    if nargin < 3
        fillval = 0;
    end
    if numel(blocksz) == 1
        blocksz = [blocksz, blocksz];
    end
    imght0 = size(IMG,1);
    imgwd0 = size(IMG,2);
    blockht = blocksz(1);
    blockwd = blocksz(2);
    Nx = ceil(imgwd0 / blockwd);
    Ny = ceil(imght0 / blockht);
    imght = blockht * Ny;
    imgwd = blockwd * Nx;
    padsz = [imght - imght0, imgwd - imgwd0];
    % IMGt = padarray(IMG, padsz, fillval, 'post');
    if any(padsz > 0)
        IMGt = fillval * ones([imght, imgwd], class(IMG));
        IMGt(1:imght0,1:imgwd0) = IMG;
    else
        IMGt = IMG;
    end
    blocks = reshape(IMGt, imght, blockwd, Nx);
    blocks = permute(blocks, [2,1,3]);
    blocks = reshape(blocks, blockwd, blockht, Nx * Ny);
    blocks = permute(blocks, [2,1,3]);
    xc = ((1 + blockwd) / 2) : blockwd : imgwd;
    yc = ((1 + blockht) / 2) : blockht : imght;
    [xx, yy] = meshgrid(xc, yc);
    cntr_yx = [yy(:), xx(:)];
end
