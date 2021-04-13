function link_struct = link_struct(links, idxs, link_struct, z_wt)
    if nargin < 4
        z_wt = [];
    end
    if nargin < 3 || isempty(link_struct)
        link_struct = [];
    end
    lnk_struct = struct;
    lnk_struct.links = links;
    lnk_struct.idx = idxs;
    lnk_struct.z_wt = z_wt;
    link_struct = [link_struct(:); lnk_struct];
end

