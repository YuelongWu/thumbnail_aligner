function s = parse_yaml(filename)
    fid = fopen(filename, 'r');
    outcell = textscan(fid, '%s%s', 'Delimiter', ':', 'Whitespace','');
    fclose(fid);
    keys = outcell{1};
    vals = outcell{2};

    cmt_marker = '#';
    keys = cellfun(@remove_comment_content, keys, ...
        repmat({cmt_marker}, size(keys)), 'UniformOutput',false);
    vals = cellfun(@remove_comment_content, vals, ...
        repmat({cmt_marker}, size(vals)), 'UniformOutput',false);
    keys = strrep(keys, char(9), repmat(' ',1,4)); % replace tab with space
    [idx, indent] = parse_indentation(keys);
    keys = strip(keys(idx));
    vals = strip(vals(idx));
    Nline = numel(keys);
    branch_flag = cellfun(@isempty, vals);
    branch_idx = find(branch_flag);
    branch_indent = indent(branch_idx);
    [branch_indent, idx] = sort(branch_indent, 'descend');
    branch_idx = branch_idx(idx);
    branch_cell = cell(numel(branch_idx), 1);
    included = false(Nline,1);
    for kb = 1:numel(branch_idx)
        b_indent = branch_indent(kb);
        b_idx = branch_idx(kb);
        t_idx = find(((1:Nline)' > b_idx) & (indent <= b_indent));
        if isempty(t_idx)
            t_idx = Nline;
        else
            t_idx = min(t_idx) - 1;
        end
        tmp = struct;
        for k = (b_idx + 1):1:t_idx
            if included(k)
                continue
            end
            if branch_flag(k)
                subtmp = branch_cell{branch_idx == k};
                tmp.(keys{k}) = subtmp;
            else
                tmp.(keys{k}) = parse_value_field(vals{k});
            end
            included(k) = true;
        end
        branch_cell{kb} = tmp;
    end
    s = struct;
    for k = 1:Nline
        if included(k)
            continue
        end
        if branch_flag(k)
            subtmp = branch_cell{branch_idx == k};
            s.(keys{k}) = subtmp;
        else
            s.(keys{k}) = parse_value_field(vals{k});
        end
    end
end

function str = remove_comment_content(str, comment_marker)
    idx = strfind(str, comment_marker);
    if ~isempty(idx)
        idx = min(idx);
        str(idx:end) = [];
    end
end

function [toKeep, indent] = parse_indentation(str)
    if ~iscell(str)
        str = {str};
    end
    N = numel(str);
    toKeep = true(N, 1);
    indent = zeros(N, 1);
    for k = 1:N
        s = str{k};
        ss = strip(s, 'left');
        if isempty(ss)
            toKeep(k) = false;
            continue
        end
        indent(k) = round((numel(s)-numel(ss))/4);
    end
    indent = indent(toKeep);
end

function val = parse_value_field(str)
    if contains(str, '[') && contains(str, ']')
        str = split(str, {',',' ','[',']',char(9),';'});
        idx = ~cellfun(@isempty, str);
        str = str(idx);
        str = replace(str, {'True','true','False','false'},{'1','1','0','0'});
        val = str2double(str).';
        return
    end
    tmp = str2double(str);
    if strcmpi(str, 'nan') || ~isnan(tmp)
        val = tmp;
    elseif strcmpi(str, 'true')
        val = true;
    elseif strcmpi(str, 'false')
        val = false;
    else
        val = str;
    end
end