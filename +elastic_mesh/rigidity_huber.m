function [costval, gradxy] = rigidity_huber(M, options)
    options = utils.set_default(options, 'huber_len', 100);
    options = utils.set_default(options, 'multiplier', 0.5);
    options = utils.set_default(options, 'min_stiff', 0);
    options = utils.set_default(options, 'compression_gain', 1);
    a = options.huber_len;
    xx = M.TR.Points(:, 1);
    yy = M.TR.Points(:, 2);
    Npt = M.pt_num;
    edg_idx = M.edges;
    dx = diff(xx(edg_idx),1,2);
    dy = diff(yy(edg_idx),1,2);
    L = (dx.^2 + dy.^2).^0.5;
    dL = L - M.L0;
    % large_mask = abs(dL) > a;
    large_mask = dL > a;
    if ~isempty(M.stiffness)
        stiffness = double(max(M.stiffness, options.min_stiff));
    else
        stiffness = 1;
    end
    costvals = 0.5 * dL.^2 .* (1 - large_mask) + ...
        a * (abs(dL) - 0.5 * a) .* large_mask;
    % penalize more on compression
    costvals(dL < 0) = options.compression_gain * costvals(dL < 0);
    costval = options.multiplier * nansum(costvals(:) .* stiffness(:));
    pgrad = dL .* (1 - large_mask) + a * sign(dL) .* large_mask;
    pgrad(dL < 0) = options.compression_gain * pgrad(dL < 0);
    L = L + 1e-5;
    pgrad_x = pgrad .* dx ./ L .* stiffness;
    pgrad_y = pgrad .* dy ./ L .* stiffness;
    grad_x = accumarray(edg_idx(:), [pgrad_x(:); -pgrad_x(:)], [Npt, 1]);
    grad_y = accumarray(edg_idx(:), [pgrad_y(:); -pgrad_y(:)], [Npt, 1]);
    gradxy = [grad_x(:), grad_y(:)] * options.multiplier;
end

