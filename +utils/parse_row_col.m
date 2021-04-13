function rc = parse_row_col(strlist, expr, rowfirst, sttpt)
    if nargin < 4
        sttpt = 1;
    end
    if nargin < 3
        rowfirst = true;
    end
    out = regexp(strlist(:),expr,'tokens');
    out = cat(1,out{:});
    out = cat(1,out{:});
    rc = str2double(out);
    if ~rowfirst
        rc = fliplr(rc);
    end
    rc = rc - sttpt + 1;
end