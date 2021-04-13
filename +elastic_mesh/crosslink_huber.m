function [costval, gradxy1, gradxy2, conf1, conf2] = crosslink_huber(links, M1, M2, options, return_grad)
    options = utils.set_default(options, 'multiplier', 1);
    options = utils.set_default(options, 'huber_len', 20);
    accumfun = [];
    if nargin < 5
        return_grad = [true, true];
    end
    a = options.huber_len;
    conf = double(links.conf);
    % xy1 = barycentricToCartesian(M1.TR,links.ID1,links.B1);
    % xy2 = barycentricToCartesian(M2.TR,links.ID2,links.B2);
    xy1 = elastic_mesh.bary2cart(M1.TR,links.ID1,links.B1);
    xy2 = elastic_mesh.bary2cart(M2.TR,links.ID2,links.B2);
    dxy = xy2 - xy1;
    dL = sum(dxy.^2, 2) .^ 0.5;
    large_mask = dL > a;
    costvals = 0.5 * dL.^2 .* (1 - large_mask) + ...
        a * (dL - 0.5 * a) .* large_mask;
    costval = options.multiplier * nansum(costvals(:) .* conf(:));
    if any(return_grad)
        pgrad = dL .* (1-large_mask) + a * large_mask;
        pgrad = pgrad .* conf;
        dL = dL + 1e-5;
        pgrad_x = pgrad .* (dxy(:,1)) ./ dL;
        pgrad_y = pgrad .* (dxy(:,2)) ./ dL;
    end
    gradxy1 = 0;
    gradxy2 = 0;
    conf1 = 1;
    conf2 = 1;
    if return_grad(1)
        cList = M1.TR.ConnectivityList;
        Npt = M1.pt_num;
        vidx = cList(links.ID1, :);
        pgrad_xb = pgrad_x .* links.B1 .* links.rmass1;
        pgrad_yb = pgrad_y .* links.B1 .* links.rmass1;
        pconf = conf .* links.B1;
        grad_x = accumarray(vidx(:), pgrad_xb(:), [Npt,1], accumfun);
        grad_y = accumarray(vidx(:), pgrad_yb(:), [Npt,1], accumfun);
        conf1 = accumarray(vidx(:), pconf(:), [Npt,1], accumfun);
        % conf1 = 1;
        gradxy1 = options.multiplier * [grad_x(:), grad_y(:)] ./ max(conf1, 0.2);
    end
    if return_grad(2)
        cList = M2.TR.ConnectivityList;
        Npt = M2.pt_num;
        vidx = cList(links.ID2, :);
        pgrad_xb = -pgrad_x .* links.B2 .* links.rmass2;
        pgrad_yb = -pgrad_y .* links.B2 .* links.rmass2;
        pconf = conf .* links.B2;
        grad_x = accumarray(vidx(:), pgrad_xb(:), [Npt,1], accumfun);
        grad_y = accumarray(vidx(:), pgrad_yb(:), [Npt,1], accumfun);
        conf2 = accumarray(vidx(:), pconf(:), [Npt,1], accumfun);
        % conf2 = 1;
        gradxy2 = options.multiplier * [grad_x(:), grad_y(:)] ./ max(conf2, 0.2);
    end
end
