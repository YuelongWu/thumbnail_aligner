function [ptpairs, IMG1, IMG2] = PMCC_matching_iter(IMG1, IMG2, maskw, maskr, options, ptpairs0, break_symmetry)
    % if either image has already been converted to descriptor, use that image as reference
    % if both converted, break the symetry
    % if none, assume the one with larger region area is the reference, biased to IMG1
    options = utils.set_default(options, 'initial_transform', 'affine');
    options = utils.set_default(options, 'initial_gridsz', 300);
    if nargin < 7
        break_symmetry = true;
    end
    if nargin < 6
        ptpairs0 = [];
    end
    cvt1 = isstruct(IMG1);
    cvt2 = isstruct(IMG2);
    switch (cvt1 + cvt2)
    case 2
        if break_symmetry
            if IMG1.ave_area >= IMG2.ave_area
                IMG2 = IMG2.img;
            else
                IMG1 = IMG1.img;
            end
            ptpairs = PMCC_matching_iter(IMG1, IMG2, maskw, maskr, options, ptpairs0);
        else
            ptpairs = PMCC_matching(IMG1, IMG2, maskw, maskr, options, ptpairs0);
        end
    case 0
        maskr1 = single(maskr{1});
        maskr2 = single(maskr{2});
        rarea1 = histcounts(maskr1(:), 0.5:1:(max(maskr1(:))+0.5));
        rarea2 = histcounts(maskr2(:), 0.5:1:(max(maskr2(:))+0.5));
        IMG10 = IMG1;
        IMG20 = IMG2;
        if mean(rarea1(rarea1 > 0)) >= mean(rarea2(rarea2 > 0))
            IMG1 = PMCC.raw_img_to_descriptor(IMG1, options, {maskw{1}, maskr{1}});
        else
            IMG2 = PMCC.raw_img_to_descriptor(IMG2, options, {maskw{2}, maskr{2}});
        end
        ptpairs = PMCC_matching_iter(IMG1, IMG2, maskw, maskr, options, ptpairs0);
        if min(sum(rarea1 > 0), sum(rarea2 > 0)) > 1
            IMG1 = IMG10;
            IMG2 = IMG20;
        end
    case 1
        if cvt1
            info1 = IMG1;
            % maskw1 = single(maskw{1});
            maskw2 = single(maskw{2});
            maskr1 = single(maskr{1});
            maskr2 = single(maskr{2});
            to_flip = false;
        else
            info1 = IMG2;
            IMG2 = IMG1;
            % maskw1 = single(maskw{2});
            maskw2 = single(maskw{1});
            maskr1 = single(maskr{2});
            maskr2 = single(maskr{1});
            ptpairs0 = geometries.flip_ptpairs(ptpairs0);
            to_flip = true;
        end

        if isempty(ptpairs0)
            maskr2 = PMCC.pad_img(maskr2, size(maskr1));
            maskrs = cat(2, maskr1(:), maskr2(:));
            idxt = all(maskrs > 0, 2);
            rgpairs = unique(maskrs(idxt,:), 'rows');
            Apm = repmat(eye(3),1,1,size(rgpairs,1));
            % Rpm = repmat(eye(3),1,1,size(rgpairs,2));
            tf_identity = true;
            options.initial_transform = 'affine';
        else
            [rgpairs, Apm, ~] = geometries.prematch2affine(ptpairs0);
            tf_identity = false;
        end
        pre_func_name = options.preprocessing.func;
        pre_func = str2func(pre_func_name);
        pre_params = options.preprocessing.params;
        IMG2f = pre_func(IMG2, pre_params, maskw2);

        if tf_identity
            IMG2f = PMCC.pad_img(IMG2f, size(maskr1));
            % maskw2 = PMCC.pad_img(maskw2, size(maskr1));
        end

        rids2 = unique(maskr2(maskr2>0));
        Nrg2 = numel(rids2);
        PTPAIRS = cell(Nrg2, 1);
        % maskr1e = mask_gen.expand_mask(maskr1);
        nb_sizes = options.detector.params.nb_size;
        Nb = numel(nb_sizes);
        if isfield(options.detector.params, 'min_dis')
            useblk = false;
            min_dises = options.detector.params.min_dis;
        else
            useblk = true;
        end

        % start aligning each sub-regions
        mask_disk = options.matcher.mask_disk;
        dis_thresh = options.matcher.dis_thresh;
        apply_mask = options.matcher.apply_mask;
        min_conf = options.matcher.min_conf;
        gridsz0 = options.matcher.gridsz;

        for kr = 1:Nrg2
            rid2 = rids2(kr);
            TFORMS = cell(numel(nb_sizes) + 1,1);
            maskr2t = (maskr2 == rid2);
            IMG2r = IMG2f .* maskr2t;
            idxr = rgpairs(:,2) == rid2;
            if all(idxr == 0)
                continue
            end
            rid1 = cast(rgpairs(idxr, 1), class(maskr1));
            maskr1t = ismember(maskr1, rid1);
            Ap = Apm(:,:,idxr);
            if strcmpi(options.initial_transform, 'affine')
                TFORMS0 = struct('rid', rgpairs(idxr, :), 'A', Ap);
            else
                idxt = ptpairs0.region_id(:,2) == rid2;
                TFORMS0 = ptpairs0(idxt, :);
            end
            TFORMS{1} = TFORMS0;
            if tf_identity
                IMG2t = IMG2r;
            else
                % IMG2t = zeros(size(maskr1), class(IMG2r));
                % for kg = 1:numel(rid1)
                %     A0 = Ap([2,1,3],[2,1,3],kg);
                %     A0 = inv(A0);
                %     A0(:,3) = [0;0;1];
                %     imgt = imwarp(IMG2r, affine2d(A0), 'linear', 'OutputView',...
                %         imref2d(size(maskr1)), 'FillValues', 0);
                %     IMG2t = IMG2t + imgt .* (maskr1 == rid1(kg));
                % end
                [~, IMG2t] = PMCC.collapse_tform_and_render(IMG2r, TFORMS, options.initial_gridsz, maskr1, rid1);
            end
            for kb = 1:Nb
                idxt = info1.kps.nb_size == nb_sizes(kb);
                if ~useblk
                    idxt = idxt & (info1.kps.min_dis == min_dises(kb));
                end
                kd = find(idxt, 1);
                kps1t = PMCC.filter_kps(info1.kps, kd, rid1);
                blk2 = PMCC.retrieve_same_block(IMG2t, info1.kps, kd, rid1);
                if dis_thresh < 1
                    dis_thresh0 = dis_thresh * nb_sizes(kb);
                else
                    dis_thresh0 = dis_thresh;
                end
                [dyx, conf] = PMCC.xcorr_fft(kps1t.des, blk2, mask_disk, [], ...
                    [kps1t.ffted, false], dis_thresh0, apply_mask);
                
                yx1 = kps1t.yx;
                yx2 = yx1 + dyx;
                region_id1 = kps1t.region_id;
                region_id = [region_id1(:), rid2 * ones(size(region_id1(:)), class(region_id1))];
                rotation = zeros(size(region_id1(:)),'single');
                ptpairst = table(yx1, yx2, region_id, conf, rotation);

                idxt = conf > min_conf;
                ptpairst = ptpairst(idxt,:);
        
                ptpairst = geometries.filter_matches_subregion(ptpairst, options.filters, nb_sizes(kb));
                TFORMS1 = ptpairst;
                if ~isempty(TFORMS1)
                    TFORMS{kb+1} = TFORMS1;
                end
                % apply transform to the image
 
                if gridsz0 > 5
                    gridsz = gridsz0;
                else
                    gridsz = gridsz0 * nb_sizes(kb);
                end
                if kb == Nb
                    if isempty(TFORMS1)
                        TFORMS0 = PMCC.collapse_tform_and_render([], TFORMS, kps1t.yx, maskr1, rid1);
                    else
                        TFORMS0 = PMCC.collapse_tform_and_render([], TFORMS, TFORMS1.yx1, maskr1, rid1);
                    end
                    TFORMS0 = TFORMS0(TFORMS0.conf > min_conf, :);
                    TFORMS0 = geometries.mask_ptpairs(TFORMS0, maskr1t, maskr2t);
                    PTPAIRS{kr} = TFORMS0;
                else
                    [~, IMG2t] = PMCC.collapse_tform_and_render(IMG2r, TFORMS, gridsz, maskr1, rid1);
                    % TFORMS0 = TFORMS0(TFORMS0.conf > min_conf, :);
                    % TFORMS0 = geometries.mask_ptpairs(TFORMS0, maskr1t, maskr2t);
                end
            end
        end
    
        idxt = cellfun(@isempty, PTPAIRS);
        PTPAIRS = PTPAIRS(~idxt);
        ptpairs = vertcat(PTPAIRS{:});
        ptpairs = geometries.filter_matches_subregion(ptpairs, options.filters, min(nb_sizes));

        if isfield(options, 'final_filters')
            ptpairs = geometries.filter_matches_subregion(ptpairs, options.final_filters, min(nb_sizes));
        end

        if to_flip
            ptpairs = geometries.flip_ptpairs(ptpairs);
            IMG2 = info1;
        end
    end
end
