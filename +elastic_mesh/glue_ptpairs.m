function ptpairs = glue_ptpairs(maskr, maskw, options)
    options = utils.set_default(options, 'wrinkle_feather', 50);
    options = utils.set_default(options, 'stiff_power', 1);
    options = utils.set_default(options, 'mesh_space', 100);

    [maskg1, maskg2] = elastic_mesh.glue_mask(maskr);

    if options.wrinkle_feather > 0
        D = bwdist(maskw);
        D = D / options.wrinkle_feather;
        maskw = min(D, 1).^(2 * options.stiff_power);
    else
        maskw = 1 - (maskw > 0);
    end

    mask_brd = maskg1 > 0;
    rg_pair = cat(2, maskg1(mask_brd), maskg2(mask_brd));
    error('need further implimentation');
end