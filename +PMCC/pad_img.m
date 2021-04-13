function IMGt = pad_img(IMG, imgsz, fillval)
    if nargin < 3
        fillval = 0;
    end
    fillval = cast(fillval, class(IMG));
    imght = size(IMG,1);
    imgwd = size(IMG,2);
    if imght < imgsz(1) && imgwd < imgsz(2)
        IMGt = IMG(1:imght, 1:imgsz, :);
    else
        idx1 = min(imgsz(1), imght);
        idx2 = min(imgsz(2), imgwd);
        imgsz0 = size(IMG);
        imgsz0(1) = imgsz(1);
        imgsz0(2) = imgsz(2);
        IMGt = fillval * ones(imgsz0, class(IMG));
        IMGt(1:idx1,1:idx2,:) = IMG(1:idx1,1:idx2,:);
    end
end
