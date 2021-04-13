function ptpairs = filter_pw_strain(ptpairs, options, varargin)
    options = utils.set_default(options, 'strain_limit', inf);
    options = utils.set_default(options, 'displc_limit', inf);
    options = utils.set_default(options, 'shear_limit', 180);
    options = utils.set_default(options, 'thresh', 0.1);
    options = utils.set_default(options, 'd_whisker', 100);
    options = utils.set_default(options, 'filter_angle', true);

    shear_limit = options.shear_limit * pi / 180;
    yx1 = ptpairs.yx1;
    yx2 = ptpairs.yx2;
    [d1, rot1] = pairwise_component_dis(yx1);
    [d2, rot2] = pairwise_component_dis(yx2);

    strain = abs(log(d1./d2));
    if isnan(options.strain_limit)
        strain_limit = 2 * median(strain(:), 'omitnan');
    else
        strain_limit = options.strain_limit;
    end
    valid_strain = strain <= strain_limit;

    displc = abs(d1 - d2);
    if isnan(options.displc_limit)
        displc_limit = median(displc(valid_strain), 'omitnan');
    else
        displc_limit = options.displc_limit;
    end
    
    valid_displc = displc <= displc_limit;
    
    valid_strain = valid_strain & valid_displc;
    
    if options.filter_angle
        rot = rot1 - rot2;
        crot = cos(rot);
        srot = sin(rot);
        crot_med = median(crot(valid_strain), 'omitnan');
        srot_med = median(srot(valid_strain), 'omitnan');
        rot_med = atan2(srot_med, crot_med);
        valid_rot = abs(wrapToPi(rot - rot_med)) <= shear_limit;
        valid_strain = valid_rot & valid_strain;
    end
    d_wt = 1./(max(d1, d2) + options.d_whisker).^2;
    valid_num = sum(valid_strain.*d_wt, 1)./sum(d_wt,1);
    valid_idx = valid_num > (options.thresh * max(valid_num));
    ptpairs = ptpairs(valid_idx, :);
end

function [d, rot] = pairwise_component_dis(yx)
    N = size(yx, 1);
    y = repmat(yx(:,1), 1, N);
    x = repmat(yx(:,2), 1, N);
    dx = x.' - x;
    dy = y.' - y;
    d = (dx.^2 + dy.^2).^0.5;
    rot = atan2(dy, dx);
end
