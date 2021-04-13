function combine_tiles(secdir, options)
    imglist = dir([secdir, filesep, '*.png']);
    imgnames = {imglist.name};
    outcell = split(imgnames(:), '_tile');
    [outname, ~, ic] = unique(outcell(:, 1));
    for k = 1 : numel(outname)
        idxt = find(ic == k);
        if numel(idxt) == 1
            movefile([imglist(idxt).folder, filesep, imglist(idxt).name], ...
                [imglist(idxt).folder, filesep, outname{k}, '.png']);
            continue
        end
        IMGt = zeros(options.block_size, options.block_size, 'single');
        WT = zeros(options.block_size, options.block_size, 'single');
        for t0 = 1 : numel(idxt)
            t = idxt(t0);
            IMG = imread([imglist(t).folder, filesep, imglist(t).name]);
            wt = local_generate_mask(IMG == options.fillval, options);
            IMGt = single(IMG) .* wt + IMGt;
            WT = wt + WT;
            delete([imglist(t).folder, filesep, imglist(t).name])
        end
        IMGout = uint8(IMGt ./ WT);
        imwrite(IMGout, [imglist(t).folder, filesep, outname{k}, '.png']);

    end
end

function wt = local_generate_mask(wt, options)
    if options.weight_power == 0
        wt = 1 - wt;
        return;
    end
    wt = bwdist(wt) + 1;
    wt = wt / options.blend_range;
    wt = min(wt,1);
    wt = wt .^ options.weight_power;
end
