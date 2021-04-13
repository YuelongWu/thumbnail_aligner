function [Ms, cost_history] = optimize_wrapper_multires(imglist, ptpairlist, fixed_sec, options, maskwdir, maskrdir, meshdir, z_wt, M0t)
    Nimg = numel(imglist);
    if Nimg == 0
        Ms = cell(0,1);
        return
    end
    if nargin < 9
        M0t = [];
    end
    if nargin < 8
        z_wt = [];
    end
    if nargin < 7
        meshdir = [];
    end
    if nargin < 6
        maskrdir = [];
    end
    if nargin < 5
        maskwdir = [];
    end
    if nargin < 4
        options = [];
    end
    if nargin < 3 || isempty(fixed_sec)
        fixed_sec = false(numel(imglist),1);
    end
    options = utils.set_default(options, 'share_M0', false);
    options = utils.set_default(options, 'rigid_first', true);
    options = utils.set_default(options, 'optimize_angle', false);
    options = utils.set_default(options, 'ext', '.png');
    options = utils.set_default(options, 'Ntranslation', 0);
    options = utils.set_default(options, 'Nrigid', 200);
    options = utils.set_default(options, 'Naffine', 200);
    options = utils.set_default(options, 'Nshared', 0);
    options = utils.set_default(options, 'Nelastic', 500);
    options = utils.set_default(options, 'lr_scl', 5);
    options = utils.set_default(options, 'milestone', 5);
    options = utils.set_default(options, 'maximum_gain', 2);
    options = utils.set_default(options, 'post_affine_filter_thresh', 0);
    h0 = options.params.learning_rate;
    % filter ptpairs list
    IMGnames = erase({imglist.name}, options.ext);

    if isfield(ptpairlist, 'imgnames')
        imgnames = erase(vertcat(ptpairlist.imgnames), options.ext);
        idxt = ismember(imgnames(:,1), IMGnames) & ismember(imgnames(:,2), IMGnames);
        ptpr_struct = ptpairlist(idxt);
    else % ptpair info is in file
        ptpr_struct = cell(numel(ptpairlist),1);
        t0 = 1;
        for k = 1:numel(ptpairlist)
            if contains(ptpairlist(k).name, '_chunk')
                load([ptpairlist(k).folder, filesep, ptpairlist(k).name], 'PTPAIRS');
                for t = 1:numel(PTPAIRS)
                    strct = PTPAIRS(t);
                    imgnames = erase(strct.imgnames, options.ext);
                    if any(strcmpi(IMGnames, imgnames{1})) && any(strcmpi(IMGnames, imgnames{2}))
                        ptpr_struct{t0} = strct;
                        t0 = t0 + 1;
                    end
                end
            else
                imgnames = utils.split_pair_name(erase({ptpairlist(k).name}, '.mat'));
                if any(strcmpi(IMGnames, imgnames{1})) && any(strcmpi(IMGnames, imgnames{2}))
                    load([ptpairlist(k).folder, filesep, ptpairlist(k).name], 'ptpairs', 'imgnames');
                    strct = struct;
                    strct.imgnames = imgnames;
                    strct.ptpairs = ptpairs;
                    ptpr_struct{t0} = strct;
                    t0 = t0 + 1;
                end
            end
        end
        idxt = ~cellfun(@isempty, ptpr_struct);
        ptpr_struct = ptpr_struct(idxt);
        ptpr_struct = vertcat(ptpr_struct{:});
    end
    MATname_split = erase(vertcat(ptpr_struct.imgnames), options.ext);
    MATnames = strcat(MATname_split(:,1), {'_'}, MATname_split(:,2));

    if ~isempty(z_wt)
        z_wt = z_wt(idxt);
    else
        z_wt = ones(numel(ptpr_struct),1, 'single');
    end
    Nlink = numel(ptpr_struct);

    % initialize mesh grid
    Ms = cell(Nimg,1);
    fixed_idx = find(fixed_sec(:));
    if options.rigid_first
        A = repmat(eye(3),1,1,Nimg);
        partnerst = cell(Nimg,1);
        matidxt = cell(Nimg,1);
    end

    for km = 1:numel(options.mesh_space)
        if options.share_M0
            refidx = max([1; find(fixed_sec(:))]);
            IMG = imread([imglist(refidx).folder, filesep, imglist(refidx).name]);
            M0 =  elastic_mesh.gen_eqtriang_mesh_mask(ones(size(IMG)), options.mesh_space(km), 1);
        end
        if km > 1
            for k = 1:Nimg
                if fixed_sec(k)
                    continue
                end
                fname = [maskwdir, filesep, imglist(k).name];
                if exist(fname, 'file')
                    maskw = imread(fname);
                else
                    maskw = [];
                end
                fname = [maskrdir, filesep, imglist(k).name];
                if exist(fname, 'file')
                    maskr = imread(fname);
                else
                    maskr = [];
                end
                if options.share_M0
                    Ms{k} = elastic_mesh.refine_mesh(Ms{k}, M0, maskr, maskw == 200, options);
                else
                    Ms{k} = elastic_mesh.refine_mesh(Ms{k}, options.mesh_space(km), maskr, maskw == 200, options);
                end
            end
        else
            for k = 1:Nimg
                if fixed_sec(k)
                    if isempty(M0t)
                        if exist([meshdir, filesep, strrep(imglist(k).name, options.ext,'.mat')], 'file')
                            load([meshdir, filesep, strrep(imglist(k).name, options.ext,'.mat')], 'M');
                            Ms{k} = M;
                        else
                            IMG = imread([imglist(k).folder, filesep, imglist(k).name]);
                            fname = [maskwdir, filesep, imglist(k).name];
                            if exist(fname, 'file')
                                maskw = imread(fname);
                            else
                                maskw = zeros(size(IMG), 'uint8');
                            end
                            fname = [maskrdir, filesep, imglist(k).name];
                            if exist(fname, 'file')
                                maskr = imread(fname);
                            else
                                maskr = ones(size(IMG), 'uint8');
                            end
                            if ~options.share_M0
                                M0 =  elastic_mesh.gen_eqtriang_mesh_mask(maskw<255, options.mesh_space(km) * 2, 1);
                            end
                            Ms{k} = elastic_mesh.init_mesh_subregion(M0, maskr, maskw == 200, options);
                        end
                    elseif numel(M0t) == 1
                        Ms{k} = M0t;
                    else
                        Ms{k} = M0t{fixed_idx == k};
                    end
                    if options.rigid_first
                        M = Ms{k};
                        xy0t = M.TR0.Points;
                        xy1t = M.TR.Points;
                        Nrgn = size(xy1t,1) / size(xy0t,1);
                        xy0t = repmat(xy0t, Nrgn, 1);
                        [~, R] = geometries.fit_affine(xy1t, xy0t);
                        A(:,:,k) = R;
                    end
                else
                    IMG = imread([imglist(k).folder, filesep, imglist(k).name]);
                    fname = [maskwdir, filesep, imglist(k).name];
                    if exist(fname, 'file')
                        maskw = imread(fname);
                    else
                        maskw = zeros(size(IMG), 'uint8');
                    end
                    fname = [maskrdir, filesep, imglist(k).name];
                    if exist(fname, 'file')
                        maskr = imread(fname);
                    else
                        maskr = ones(size(IMG), 'uint8');
                    end
                    if ~options.share_M0
                        M0 =  elastic_mesh.gen_eqtriang_mesh_mask(maskw<255, options.mesh_space(km), 1);
                    end
                    Ms{k} = elastic_mesh.init_mesh_subregion(M0, maskr, maskw == 200, options);
                end

                if options.rigid_first
                    imgnm = strrep(imglist(k).name, options.ext, '');
                    idxt = find(ismember(MATname_split(:,1), imgnm) | ismember(MATname_split(:,2),imgnm));
                    matnames = MATnames(idxt);
                    matnames = strrep(matnames, imgnm, '');
                    [~, k1s] = ismember(strrep(matnames, '_', ''), strrep(IMGnames,'_',''));
                    matidxt{k} = idxt;
                    partnerst{k} = k1s;
                end
            end

            if options.rigid_first
                psz = cellfun(@numel, partnerst);
                partners = zeros(Nimg, max(psz(:)), 'uint32');
                matidx = zeros(Nimg, max(psz(:)), 'uint32');
                for k = 1:Nimg
                    k1s = partnerst{k};
                    partners(k, 1:numel(k1s)) = k1s;
                    matidx(k, 1:numel(k1s)) = matidxt{k};
                end
                processed_nb = ismember(partners, fixed_idx);
                to_align = ~fixed_sec(:);
                while any(to_align)
                    resv_ref = sum(processed_nb, 2) .* to_align;
                    if all(resv_ref == 0)
                        % no connected pairs remaining, start with one with most neighbor
                        [~, k0] = max((psz(:)+1) .* to_align(:));
                    else
                        [n_ref, k0] = max(resv_ref);
                        idxt = matidx(k0, :) .* cast(processed_nb(k0,:), class(matidx));
                        idxt = idxt(idxt > 0);
                        XY0 = cell(n_ref,1);
                        XY1 = cell(n_ref,1);
                        for t = 1:n_ref
                            ptpairs = ptpr_struct(idxt(t)).ptpairs;
                            imgnames = erase(ptpr_struct(idxt(t)).imgnames, options.ext);
                            k1 = find(strcmpi(IMGnames, imgnames{1}));
                            k2 = find(strcmpi(IMGnames, imgnames{2}));
                            if k1 == k0
                                ptpairs = geometries.flip_ptpairs(ptpairs);
                                kt = k2;
                            else
                                kt = k1;
                            end
                            XY1{t} = fliplr(ptpairs.yx2);
                            XY0{t} = fliplr(ptpairs.yx1) * A(1:2,1:2,kt) + A(3,1:2,kt);
                        end
                        xy1 = cat(1, XY1{:});
                        xy0 = cat(1, XY0{:});
                        [~, R] = geometries.fit_affine(xy0, xy1);
                        A(:,:,k0) = R;
                    end
                    to_align(k0) = false;
                    processed_nb(partners == k0) = 1;
                end
                for t = 1:Nimg
                    if ~fixed_sec(t)
                        R = A(:,:,t);
                        Ms{t}.TR.Points = Ms{t}.TR.Points * R(1:2,1:2) + R(3,1:2);
                    end
                end
            end
        end

        links =  cell(Nlink,1);
        for k = 1:Nlink
            ptpairs = ptpr_struct(k).ptpairs;
            imgnames = ptpr_struct(k).imgnames;
            k1 = find(strcmpi(IMGnames, erase(imgnames{1}, options.ext)));
            k2 = find(strcmpi(IMGnames, erase(imgnames{2}, options.ext)));
            lnk = elastic_mesh.ptpairs_to_links(ptpairs, Ms{k1}, Ms{k2});
            lnk_struct = elastic_mesh.link_struct(lnk, [k1, k2], [], z_wt(k));
            links{k} = lnk_struct;
        end
        links = vertcat(links{:});

        Naffine = options.Naffine;
        Nelastic = options.Nelastic;
        Nshared = options.Nshared;
        lr_scl = options.lr_scl;
        cost_history = cell(Naffine + Nelastic + Nshared + options.Nrigid, 1);
        options.params.learning_rate = h0 * lr_scl;
        t = 1;
        for k = 1:options.Ntranslation
            [Ms, cost_history0] = elastic_mesh.BGD(Ms, links, options.params, fixed_sec, 'GLOBAL_TRANSLATION');
        end
        if options.Ntranslation > 0 && options.post_translate_filter_thresh > 0
            DD = cell(Nlink, 1);
            dmax = 0;
            for k = 1:Nlink
                lnk_struct = links(k);
                sec_indices = lnk_struct.idx;
                lnk = lnk_struct.links;
                M1 = Ms{sec_indices(1)};
                M2 = Ms{sec_indices(2)};
                xy1t = elastic_mesh.bary2cart(M1.TR, lnk.ID1, lnk.B1);
                xy2t = elastic_mesh.bary2cart(M2.TR, lnk.ID2, lnk.B2);
                dis = sum((xy1t - xy2t) .^ 2, 2) .^ 0.5;
                DD{k} = dis;
                dmax = max(dmax, max(dis(:)));
            end
            to_keep = true(Nlink,1);
            if dmax > options.post_translate_filter_thresh
                d_thresh = max(dmax/2, options.post_translate_filter_thresh);
                for k = 1:Nlink
                    discrd = DD{k} > d_thresh;
                    if any(discrd(:))
                        if all(discrd(:))
                            to_keep(k) = false;
                        else
                            lnk_struct = links(k);
                            lnk_struct.links = lnk_struct.links(~discrd, :);
                            links(k) = lnk_struct;
                        end
                    end
                end
                links = links(to_keep);
                Nlink = numel(links);
            end
        end
        for k = 1:options.Nrigid
            [Ms, cost_history0] = elastic_mesh.BGD(Ms, links, options.params, fixed_sec, 'REGION_RIGID');
        end
        cost_history{t} = cost_history0(~isnan(cost_history0));

        if Nshared > 0
            options.params.learning_rate = h0;
            h1 = options.params.learning_rate;
            c1 = nan;
            cnt = 0;
            cnt_suc = 0;
            for k = 1:Nshared
                t = t+1;
                [Ms, cost_history0] = elastic_mesh.BGD(Ms, links, options.params,fixed_sec, 'SHARED_ELASTIC');
                cost_history{t} = cost_history0;
                if k == 1
                    c1 = cost_history0(end);
                else
                    c0 = c1;
                    c1 = cost_history0(end);
                    if c1 >= c0
                        cnt = cnt + 1;
                        cnt_suc = 0;
                    else
                        cnt = 0;
                        cnt_suc = cnt_suc + 1;
                    end
                end
                if cnt > options.milestone
                    cnt = cnt - 1;
                    options.params.learning_rate = options.params.learning_rate / 1.25;
                    if options.params.learning_rate < 0.25 * h1
                        break
                    end
                end
                if cnt_suc > options.milestone
                    options.params.learning_rate = options.params.learning_rate * 1.25;
                    if options.params.learning_rate > h1 * options.maximum_gain
                        options.params.learning_rate = h1 * options.maximum_gain;
                    end
                end
            end
            if options.params.learning_rate > h1
                warning([imglist(1).name, '->',imglist(end).name,':Shared elastic optimization stopped early..']);
            end
            % relax rigid spring
            Ms = cellfun(@elastic_mesh.relax_rigidity_spring, Ms(:), 'UniformOutput', false);
        end

        if Naffine > 0
            options.params.learning_rate = h0 * lr_scl;
            h1 = options.params.learning_rate;
            c1 = nan;
            cnt = 0;
            cnt_suc = 0;
            for k = 1:Naffine
                t = t+1;
                % tic;
                [Ms, cost_history0] = elastic_mesh.BGD(Ms, links, options.params,fixed_sec, 'REGION_AFFINE');
                cost_history{t} = cost_history0;
                if k == 1
                    c1 = cost_history0(end);
                    % toc;
                else
                    c0 = c1;
                    c1 = cost_history0(end);
                    if c1 >= c0
                        cnt = cnt + 1;
                        cnt_suc = 0;
                    else
                        cnt = 0;
                        cnt_suc = cnt_suc + 1;
                        % toc;
                    end
                end
                if cnt > options.milestone
                    options.params.learning_rate = options.params.learning_rate / 1.25;
                    if options.params.learning_rate < 0.25 * h1
                        % disp([imglist(1).name, '->',imglist(end).name,':Affine optimization converged!']);
                        break
                    end
                end
                if cnt_suc > options.milestone
                    options.params.learning_rate = options.params.learning_rate * 1.25;
                    if options.params.learning_rate > h1 * options.maximum_gain
                        options.params.learning_rate = h1 * options.maximum_gain;
                    end
                end
            end
            if options.params.learning_rate > h1
                warning([imglist(1).name, '->',imglist(end).name,':Affine optimization stopped early..']);
            end
        end

        if options.post_affine_filter_thresh > 0
            DD = cell(Nlink, 1);
            dmax = 0;
            for k = 1:Nlink
                lnk_struct = links(k);
                sec_indices = lnk_struct.idx;
                lnk = lnk_struct.links;
                M1 = Ms{sec_indices(1)};
                M2 = Ms{sec_indices(2)};
                xy1t = elastic_mesh.bary2cart(M1.TR, lnk.ID1, lnk.B1);
                xy2t = elastic_mesh.bary2cart(M2.TR, lnk.ID2, lnk.B2);
                dis = sum((xy1t - xy2t) .^ 2, 2) .^ 0.5;
                DD{k} = dis;
                dmax = max(dmax, max(dis(:)));
            end
            to_keep = true(Nlink,1);
            if dmax > options.post_affine_filter_thresh
                d_thresh = max(dmax/2, options.post_affine_filter_thresh);
                for k = 1:Nlink
                    discrd = DD{k} > d_thresh;
                    if any(discrd(:))
                        if all(discrd(:))
                            to_keep(k) = false;
                        else
                            lnk_struct = links(k);
                            lnk_struct.links = lnk_struct.links(~discrd, :);
                            links(k) = lnk_struct;
                        end
                    end
                end
                links = links(to_keep);
                Nlink = numel(links);
            end
        end

        if Nelastic > 0
            options.params.learning_rate = h0;
            h1 = options.params.learning_rate;
            c1 = nan;
            cnt = 0;
            cnt_suc = 0;
            for k = 1:Nelastic
                t = t+1;
                % tic;
                [Ms, cost_history0] = elastic_mesh.BGD(Ms, links, options.params, fixed_sec, 'ELASTIC');
                cost_history{t} = cost_history0;
                if k == 1
                        c1 = cost_history0(end);
                else
                    c0 = c1;
                    c1 = cost_history0(end);
                    if c1 >= c0
                        cnt = cnt + 1;
                        cnt_suc = 0;
                    else
                        cnt = 0;
                        cnt_suc = cnt_suc + 1;
                    end
                end
                if cnt > options.milestone
                    options.params.learning_rate = options.params.learning_rate / 1.25;
                    if options.params.learning_rate < 0.25 * h1
                        % disp([imglist(1).name, '->',imglist(end).name,':Elastic optimization converged!']);
                        break
                    end
                end
                if cnt_suc > options.milestone
                    options.params.learning_rate = options.params.learning_rate * 1.25;
                    if options.params.learning_rate > h1 * options.maximum_gain
                        options.params.learning_rate = h1 * options.maximum_gain;
                    end
                end
            end
            if options.params.learning_rate > h1
                warning([erase(imglist(1).name,options.ext), '->',erase(imglist(end).name, options.ext),':Elastic optimization stopped early..']);
            end
        end
    end

    if options.optimize_angle && ~any(fixed_sec)
        Ms = elastic_mesh.optimize_mesh_rotation(Ms);
    end
end
