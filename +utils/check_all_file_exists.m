function TF = check_all_file_exists(fpath, fname, flag)
    if nargin < 3
        flag = 'file';
    end
    fdir = strcat(fpath, filesep, fname);
    TF = true;
    for k = 1:numel(fdir)
        if ~exist(fdir{k}, flag)
            TF = false;
            return
        end
    end
end