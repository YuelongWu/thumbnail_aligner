function options = set_default(options, fieldname, defval)
    if isempty(options)
        options = struct;
    end
    if ~isfield(options, fieldname)
        options.(fieldname) = defval;
    end
end