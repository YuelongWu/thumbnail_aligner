function [Ms_out, cost_history] = BGD(Ms, links, options, fixed_mask, rigid, clip)
    % batch gradient descent
    options = utils.set_default(options, 'print_time', true);
    options = utils.set_default(options, 'max_iter', 5000);
    options = utils.set_default(options, 'learning_rate', 0.1);
    options = utils.set_default(options, 'min_learn_rate', 0);
    options = utils.set_default(options, 'momentum', 0.5);
    options = utils.set_default(options, 'momentum_diminish', 0.5);
    options = utils.set_default(options, 'smooth_gradient', false);
    options = utils.set_default(options, 'smooth_start_knsz', 1);
    options = utils.set_default(options, 'smooth_knsz_diminish', 0.75);
    options = utils.set_default(options, 'smooth_end_kernel_sz', 0.1);
    options = utils.set_default(options, 'distance_decay_power', 0);
    tic;
    Nsec = numel(Ms);
    if nargin < 6
        clip = false;
    end
    if nargin < 5
        rigid = 0;
    end
    rigid_flag_dict = {'ELASTIC'; 'GLOBAL_AFFINE'; 'GLOBAL_RIGID';...
        'REGION_AFFINE'; 'REGION_RIGID'; 'SHARED_ELASTIC'; 'GLOBAL_TRANSLATION'};
    if ~isnumeric(rigid)
        % 0 - ELASTIC (default)
        % 1 - GLOBAL_AFFINE
        % 2 - GLOBAL_RIGID
        % 3 - REGION_AFFINE
        % 4 - REGION_RIGID
        % 5 - SHARED_ELASTIC
        % 6 - GLOBAL_TRANSLATION
        rigid = find(strcmpi(rigid_flag_dict, rigid));
        if isempty(rigid)
            rigid = 0;
        else
            rigid = rigid - 1;
        end
    end
    if nargin < 4 || isempty(fixed_mask)
        fixed_mask = false(Nsec, 1);
    end
    
    Nlink = numel(links);
    h0 = options.learning_rate;
    h = h0;
    smooth_kern = Ms{1}.smooth_kernel;
    if options.smooth_gradient && ~isempty(smooth_kern)
        kensz = options.smooth_start_knsz;
        smooth_kern = smooth_kern.^(1/kensz);
    end
    mmnt = options.momentum;
    mntdim = options.momentum_diminish;
    mmnt_func = @(g1, g2) g1 + mmnt * g2;
    mmnt_func2 = @(g) mntdim * g;
    cost0 = inf;
    grad_mmnt = repmat({0}, Nsec, 1);
    Ms_out = Ms;
    cost_history = nan(options.max_iter, 1);
    sideL = inf;
    if clip
        for k = 1:numel(Ms)
            if Ms{k}.sideL < sideL
                sideL = Ms{k}.sideL;
            end
        end
    end
    for k0 = 1:options.max_iter
        cost1 = 0;
        grad_crnt = repmat({0}, Nsec, 1);
        rigid_func = str2func(options.intra_section.func);
        rigid_params = options.intra_section.params;
        for k1 = 1:Nsec
            if fixed_mask(k1)
                continue
            end
            M = Ms{k1};
            [costval, gradxy] = rigid_func(M, rigid_params);
            cost1 = cost1 + costval;
            grad_crnt{k1} = grad_crnt{k1} + gradxy;
        end
        crosslink_func = str2func(options.cross_section.func);
        clnk_params = options.cross_section.params;
        for k1 = 1:Nlink
            lnk_struct = links(k1);
            sec_indices = lnk_struct.idx;
            lnk = lnk_struct.links;
            M1 = Ms{sec_indices(1)};
            M2 = Ms{sec_indices(2)};
            compute_grad = ~fixed_mask(sec_indices);
            [costval, gradxy1, gradxy2, conf1, conf2] = crosslink_func(lnk, M1, M2, clnk_params, compute_grad);
            cost1 = cost1 + costval;
            if options.smooth_gradient && kensz >= options.smooth_end_kernel_sz
                gradxy1 = local_combine_grad_conf(gradxy1, conf1, smooth_kern);
                gradxy2 = local_combine_grad_conf(gradxy2, conf2, smooth_kern);
            end
            if options.distance_decay_power ~= 0
                dis_coeff = 1/(abs(diff(sec_indices))+1).^options.distance_decay_power;
            else
                dis_coeff = 1;
            end
            if isfield(lnk_struct, 'z_wt') && ~isempty(lnk_struct.z_wt)
                dis_coeff = dis_coeff * lnk_struct.z_wt;
            end
            grad_crnt{sec_indices(1)} = grad_crnt{sec_indices(1)} + dis_coeff * gradxy1;
            grad_crnt{sec_indices(2)} = grad_crnt{sec_indices(2)} + dis_coeff * gradxy2;
        end
        if options.smooth_gradient
            smooth_kern = smooth_kern .^ (1/options.smooth_knsz_diminish);
            kensz = kensz * options.smooth_knsz_diminish;
        end

        cost_history(k0) = cost1;
        if cost1 < cost0
            cost0 = cost1;
            Ms_out = Ms;
            grad_mmnt = cellfun(mmnt_func, grad_crnt, grad_mmnt, 'UniformOutput', false);
            if clip
                limgrad = repmat({sideL}, Nsec, 1);
                grad_mmnt = cellfun(@local_limit_grad_amp, grad_mmnt, limgrad, 'UniformOutput', false);
            end
            hc = repmat({h}, Nsec, 1);
            if rigid == 2
                Ms = cellfun(@local_apply_gradient_rigid, Ms(:), grad_mmnt(:), hc, 'UniformOutput', false);
            elseif rigid == 1
                Ms = cellfun(@local_apply_gradient_affine, Ms(:), grad_mmnt(:), hc, 'UniformOutput', false);
            elseif rigid == 4
                Ms = cellfun(@local_apply_gradient_rigid_subregion, Ms(:), grad_mmnt(:), hc, 'UniformOutput', false);
            elseif rigid == 3
                Ms = cellfun(@local_apply_gradient_affine_subregion, Ms(:), grad_mmnt(:), hc, 'UniformOutput', false);
            elseif rigid == 6
                Ms = cellfun(@local_apply_gradient_translation, Ms(:), grad_mmnt(:), hc, 'UniformOutput', false);
            else
                if rigid == 5
                    grad_mmnt0 = 0;
                    grad_cnt0 = 0;
                    for k1 = 1:Nsec
                        if fixed_mask(k1)
                            continue
                        end
                        grad_mmnt0 = grad_mmnt0 + grad_mmnt{k1};
                        grad_cnt0 = grad_cnt0 + 1;
                    end
                    grad_mmnt0 = grad_mmnt0 / grad_cnt0;
                    grad_mmnt(~fixed_mask) = {grad_mmnt0}; 
                end
                Ms = cellfun(@local_apply_gradient, Ms(:), grad_mmnt(:), hc, 'UniformOutput', false);
            end
            h = h * 1.1;
            % M = Ms_out{find(~fixed_mask, 1)};xy0 = repmat(M.TR0.Points,M.region_num);xy1 = M.TR.Points;quiver(xy0(:,1), xy0(:,2), xy1(:,1)-xy0(:,1),xy1(:,2)-xy0(:,2));
            % M = Ms_out{find(~fixed_mask, 1)};xy0 = repmat(M.TR0.Points,M.region_num);xy1 = M.TR.Points;A = geometries.fit_affine(xy0, xy1); xy1=xy1*A(1:2,1:2)+A(3,1:2);quiver(xy0(:,1), xy0(:,2), xy1(:,1)-xy0(:,1),xy1(:,2)-xy0(:,2),'autoscale','off');
            if h > 6 * h0
                h = 6 * h0;
            end
        else
            Ms = Ms_out;
            h = h * 0.5;
            % grad_mmnt = repmat({0}, Nsec, 1);
            grad_mmnt = cellfun(mmnt_func2, grad_mmnt, 'UniformOutput', false);
        end
        if h < options.min_learn_rate * options.learning_rate
            break
        end
    end
    cost_history = cost_history(~isnan(cost_history));
    if options.print_time
        t = toc;
        disp([rigid_flag_dict{rigid + 1}, ': ', num2str(t),' seconds']);
    end
end

function gradxy = local_limit_grad_amp(gradxy, thresh)
    d = sum(gradxy.^2,2).^0.5;
    idx = d > thresh;
    if any(idx)
        gradxy(idx, :) = gradxy(idx, :)./ d(idx) * thresh;
    end
end

function gradxy = local_combine_grad_conf(gradxy, conf, smooth_kern)
    % smooth the gradient
    if isempty(smooth_kern)
        return
    end
    conf = double(conf(:));
    gradxy = gradxy .* conf;
    Npt = size(gradxy, 1);
    Npt0 = size(smooth_kern, 1);
    Nrg = floor(Npt/Npt0);
    t0 = 0;
    for k = 1:Nrg
        gradxy((t0+1):(t0+Npt0),:) = smooth_kern * gradxy((t0+1):(t0+Npt0),:);
        conf((t0+1):(t0+Npt0)) = smooth_kern * conf((t0+1):(t0+Npt0));
    end
    conf = max(conf, 1e-6);
    gradxy = gradxy./conf;
end

function M = local_apply_gradient(M, grad, h)
    M.TR.Points = M.TR.Points + h * grad;
end

function M = local_apply_gradient_translation(M, grad, h)
    idxt = any(abs(grad) > 0, 2);
    if any(idxt)
        mxy = mean(grad(idxt, :), 1);
        M.TR.Points = M.TR.Points + h * mxy;
    end
end


function M = local_apply_gradient_rigid(M, grad, h)
    tmp = M.TR.Points + h * grad;
    [~, R] = geometries.fit_affine(tmp, M.TR.Points);
    % R = geometries.fit_affine(tmp, M.TR.Points);
    M.TR.Points = M.TR.Points * R(1:2,1:2) + R(3,1:2);
end

function M = local_apply_gradient_affine(M, grad, h)
    tmp = M.TR.Points + h * grad;
    R = geometries.fit_affine(tmp, M.TR.Points);
    M.TR.Points = M.TR.Points * R(1:2,1:2) + R(3,1:2);
end

function M = local_apply_gradient_rigid_subregion(M, grad, h)
    Nrg = M.region_num;
    pt_num = M.pt_num0;
    if all(grad(:) == 0)
        return
    end
    tmp = M.TR.Points + h * grad;
    for g = 1:Nrg
        idx0 = 1 + (g - 1) * pt_num;
        idx1 = g * pt_num;
        [~,R] = geometries.fit_affine(tmp(idx0:idx1,:), M.TR.Points(idx0:idx1,:));
        M.TR.Points(idx0:idx1,:) = M.TR.Points(idx0:idx1,:) * R(1:2,1:2) + R(3,1:2);
    end
end

function M = local_apply_gradient_affine_subregion(M, grad, h)
    Nrg = M.region_num;
    pt_num = M.pt_num0;
    if all(grad(:) == 0)
        return
    end
    tmp = M.TR.Points + h * grad;
    for g = 1:Nrg
        idx0 = 1 + (g - 1) * pt_num;
        idx1 = g * pt_num;
        R = geometries.fit_affine(tmp(idx0:idx1,:), M.TR.Points(idx0:idx1,:));
        M.TR.Points(idx0:idx1,:) = M.TR.Points(idx0:idx1,:) * R(1:2,1:2) + R(3,1:2);
    end
end