function arrayout = padones(arrayin, dim)
    if nargin < 2
        dim = 2;
    end
    sz = size(arrayin);
    sz(dim) = 1;
    addone = ones(sz);
    arrayout = cat(dim, arrayin, addone);
end