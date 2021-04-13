function filters = init_radon_filt(options)
    filters = struct;
    filters = copy_settings(options, filters, 'proj_num', 12);
    filters = copy_settings(options, filters, 'beam_num', 15);
    filters = copy_settings(options, filters, 'beam_wd', 5);
    filters = copy_settings(options, filters, 'eccentric', false);
    filters = copy_settings(options, filters, 'filter_blur', 0);
    filters = copy_settings(options, filters, 'beam_radius', ceil(options.beam_num * options.beam_wd / 2));

    proj_num = filters.proj_num;
    beam_num = filters.beam_num;
    angles = ((1:proj_num) - 1) * 180 / proj_num;
    filters.kernels = cell(proj_num, 1);
    offset0 = zeros(beam_num, 2);
    offset0(:, 2) = filters.beam_wd * linspace(1 - beam_num, beam_num - 1, beam_num) / 2;
    if filters.eccentric
        kernel0 = ones(filters.beam_radius, 1);
        offset0(:, 1) = filters.beam_radius/2;
    else
        kernel0 = ones(2*filters.beam_radius+1, 1);
    end
    offsets = repmat(offset0, 1, 1, 2 * proj_num);

    for k = 1:proj_num
        kernel1 = imrotate(kernel0, angles(k), 'bilinear');
        filters.kernels{k} = kernel1 / sum(kernel1(:));
        cs = cosd(angles(k));
        sn = sind(angles(k));
        R = [cs, sn; -sn, cs];
        offset1 = offset0 * R;
        offsets(:,:,k) = offset1;
        offsets(:,:,k+proj_num) = -offset1;
    end
    filters.offsets = offsets;
    filters.beam_se = beam_thicken_disk(options.beam_wd);
end

function filters = copy_settings(options, filters, fieldname, defval)
    if isfield(filters, fieldname)
        return
    end
    if ~isfield(options, fieldname)
        filters.(fieldname) = defval;
    else
        filters.(fieldname) = options.(fieldname);
    end
end

function se = beam_thicken_disk(beam_wd)
    if beam_wd < 2
        se = [];
        return
    end
    sideL = 2 * floor(beam_wd/2) + 1;
    mask = false(sideL, sideL);
    mask((numel(mask)+1)/2) = true;
    D = bwdist(mask);
    d = D(:);
    r1 = (1/pi)^0.5;
    r2 = beam_wd/2;
    d1 = (r1^2 - r2^2 + d.^2)./(2*d);
    d2 = d - d1;
    dA = r1^2*acos(d1/r1) - d1.*(r1^2-d1.^2).^0.5 + r2^2*acos(d2/r2) - d2.*(r2^2 - d2.^2).^0.5;
    dA = real(dA)/(pi*r1^2);
    dA = max(min(dA, 1),0);
    se = reshape(dA, size(D));
    se = double(se / nansum(se(:)));
end