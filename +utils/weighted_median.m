function m = weighted_median(A, WT, dim, q)
    if nargin < 4
        q = 0.5;
    end
    szs = size(A);
    Ndim = max(dim, numel(szs));
    pmt = unique([dim, 1:Ndim], 'stable');
    A = permute(A, pmt);
    WT = permute(WT, pmt);
    szsn = size(A);
    Ax = single(A) + single(WT) * 1i;
    Ax = sort(Ax, 1, 'ComparisonMethod', 'real');
    stcksz = size(A, 1);
    mWT = min(WT(WT > 0));
    mWT = mWT / (100 * stcksz);
    WTx = max(mWT, imag(Ax));
    WTx = WTx ./ sum(WTx, 1);
    WTsum = cumsum(WTx, 1);
    idxt = ((WTsum - q) .* (WTsum - WTx - q)) <= 0;
    A = real(Ax) .* idxt;
    m = sum(A, 1) ./ sum(idxt, 1);
    m = reshape(m, [szsn(2:end),1]);
end
