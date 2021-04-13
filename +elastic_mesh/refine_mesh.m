function Mt = refine_mesh(M, mesh_space, maskr, maskw, options)
    if nargin < 5
        options = struct;
    end
    if ~isstruct(mesh_space)
        if M.sideL == mesh_space
            Mt = M;
            return
        end

        outsz = mesh_space * ceil(max(fliplr(M.TR0.Points)) / mesh_space + 1);
        if nargin < 4 || isempty(maskw)
            maskw = [];
        else
            outsz = size(maskw);
        end
        if nargin < 3 || isempty(maskr)
            maskr = [];
        else
            outsz = size(maskr);
        end
        if isempty(maskw)
            maskw = zeros(outsz);
        end
        if isempty(maskr)
            yx0 = elastic_mesh.gen_eqtriang_mesh(outsz, mesh_space, false, false);
            yx0 = max(min(round(yx0), outsz), 1);
            maskr = zeros(outsz);
            ind0 = sub2ind(outsz, yx0(:,1), yx0(:,2));
            maskr(ind0) = 1;
            [~, indx] = bwdist(maskr);
            ID = pointLocation(M.TR0, fliplr(yx0));
            maskr(ind0(isnan(ID))) = 0;
            maskr = maskr(indx);
        end
        M0 = elastic_mesh.gen_eqtriang_mesh_mask(maskr > 0, mesh_space, 1);
    else
        M0 = mesh_space;
        if M.sideL == M0.sideL
            Mt = M;
            return
        end 
    end
    Mt = elastic_mesh.init_mesh_subregion(M0, maskr, maskw, options);
    [B, ID] = elastic_mesh.cart2bary(M.TR0, fliplr(Mt.TR0.Points));
    for k = 1 : Mt.region_num
        idx1 = 1 + (k-1) * Mt.pt_num0;
        idx2 = k * Mt.pt_num0;
        offst = (k-1) * M.pt_num0;
        Mt.TR.Points(idx1:idx2,:) = elastic_mesh.bary2cart(M.TR, ID, B, offst);
    end
end
