function parsave(filename, x, varname)
    if iscell(varname)
        for k = 1:numel(varname)
            t = x{k}; %#ok
            eval([varname{k}, ' = t;']);
        end
        save(filename, varname{:});
    else
        eval([varname, ' = x;']);
        save(filename, varname);
    end
end