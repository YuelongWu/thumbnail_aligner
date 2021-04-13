function x = parload(fname, varname)
    load(fname, varname);
    eval(['x = ', varname, ';']);
end
