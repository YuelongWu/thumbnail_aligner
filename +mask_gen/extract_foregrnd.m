function [mask, thresh] = extract_foregrnd(IMG, mask, bordersz, qt)
    if nargin < 4
        qt = 0.75;
    end
    if nargin < 3 || isempty(bordersz)
        bordersz = [5, 15];
    end
    if numel(bordersz) == 1
        bordersz = [0, bordersz];
    end
    if nargin < 2 || isempty(mask)
        mask = IMG == 255;
    end
    if isinf(bordersz(2))
        idxt1 = true(size(mask));
    else
        idxt1 = imdilate(mask, strel('disk', bordersz(2)));
    end
    if bordersz(1) == 0
        idxt2 = ~mask;
    else
        idxt2 = ~imdilate(mask, strel('disk', bordersz(1)));
    end 
    idxt = idxt1 & idxt2;
    thresh = single(quantile(IMG(idxt), qt));
    mask = IMG < (2*thresh-255);
end
