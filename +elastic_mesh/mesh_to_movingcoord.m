function D = mesh_to_movingcoord(TR0, TR1, imgsz, meshsize)
    if nargin < 4
        meshsize = [100, 100];
    end
    if numel(meshsize) == 1
        meshsize = [meshsize, meshsize];
    end
    imght = imgsz(1);
    imgwd = imgsz(2);
    Nx = ceil(imgwd / meshsize(2));
    Ny = ceil(imght / meshsize(1));
    gx = (0:Nx) * meshsize(2);
    gx = gx - mean(gx(:)) + (1 + imgwd) / 2;
    gy = (0:Ny) * meshsize(1);
    gy = gy - mean(gy(:)) + (1 + imght) / 2;
    [gxx, gyy] = meshgrid(gx, gy);
    [B, ID] = elastic_mesh.cart2bary(TR1, [gyy(:), gxx(:)]);
    gxy1 = barycentricToCartesian(TR0,ID,B);
    gxx1 = gxy1(:,1);
    gyy1 = gxy1(:,2);
    DX = interp2(gxx, gyy, reshape(gxx1, size(gxx)),...
        (1:imgwd)',1:imght, 'makima'); % makima
    DY = interp2(gxx, gyy, reshape(gyy1, size(gyy)),...
        (1:imgwd)',1:imght, 'makima');
    D = cat(3, DX, DY);
end