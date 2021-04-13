function [idx, isseed] = filter_img_from_mat(IMGname_all, MATname, IMGname_seed, ext)
    if nargin < 4
        ext = '.png';
    end
    MATname = erase(MATname, '.mat');
    IMGname_all = erase(IMGname_all, ext);
    if nargin < 3 || isempty(IMGname_seed)
        MATname_split = utils.split_pair_name(MATname);
        idx = ismember(IMGname_all(:), MATname_split(:));
        idx = find(idx);
        isseed = false(size(idx));
        return
    end
    MATname_split = utils.split_pair_name(MATname);
    IMGname_seed = erase(IMGname_seed, ext);
    idxt = ismember(MATname_split(:,1),IMGname_seed) | ismember(MATname_split(:,2),IMGname_seed);
    MATname = MATname(idxt);
    MATname_split = utils.split_pair_name(MATname);
    idx = ismember(IMGname_all(:), MATname_split(:));
    idx = find(idx);
    isseed = ismember(IMGname_all(idx), IMGname_seed);
end
