function ptpairs = filter_matches(ptpairs, options, L)
    if nargin < 3
        L = [];
    end
    if isfield(options, 'filters')
        options = options.filters;
    end
    filternames = fieldnames(options);
    filteridx = contains(filternames, 'filter');
    filternames = filternames(filteridx);
    for k = 1:numel(filternames)
        fltr = options.(filternames{k});
        fltr_name = fltr.func;
        fltr_func = str2func(fltr_name);
        fltr_params = fltr.params;
        ptpairs = fltr_func(ptpairs, fltr_params, L);
    end
end