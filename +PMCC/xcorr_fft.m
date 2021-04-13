function [dyx , conf] = xcorr_fft(IMG1, IMG2, mask_disk, expected_dyx, FFTed, dis_thresh, apply_mask, xcorr_normalize)
    if nargin < 8
        xcorr_normalize = 0;
    end
    if nargin < 7
        apply_mask = false;
    end
    if nargin < 6
        dis_thresh = inf;
    end
    if dis_thresh < 0
        flipcmp = true;
        dis_thresh = abs(dis_thresh);
    else
        flipcmp = false;
    end
    if dis_thresh < 1
        dis_thresh = dis_thresh * min(size(IMG1,1), size(IMG1,2));
    end
    if nargin < 5
        FFTed = [false, false];
    end
    if isscalar(FFTed)
        FFTed = [FFTed, FFTed];
    end
    if mask_disk < 1
        mask_disk = mask_disk * max(size(IMG1,1), size(IMG1,2));
        mask_disk = max(mask_disk, 4);
    end
    szshift = ([size(IMG2, 1), size(IMG2,2)] - [size(IMG1, 1), size(IMG1,2)]) / 2;
    imght = max(size(IMG1, 1), size(IMG2, 1));
    imgwd = max(size(IMG1, 2), size(IMG2, 2));
    Nimg = max(size(IMG1, 3), size(IMG2, 3));
    if nargin < 4 || isempty(expected_dyx)
        expected_dyx = zeros(Nimg, 2);
    end
    expected_dyx = expected_dyx + szshift;
    if ~FFTed(1)
        IMG1 = fft2(IMG1, imght, imgwd);
    end
    if ~FFTed(2)
        IMG2 = fft2(IMG2, imght, imgwd);
    end
    F = conj(IMG1) .* IMG2;
    C0 = real(ifft2(F));
    if xcorr_normalize > 0
        msk = false(size(C0));
        if xcorr_normalize == 1
            msk([1,end],:) = 1;
            D = bwdist(msk).^0.5;
        elseif xcorr_normalize == 2
            msk = false(size(C0));
            msk(:, [1,end]) = 1;
            D = bwdist(msk).^0.5;
        else
            msk = false(size(C0));
            msk(:, [1,end]) = 1;
            D1 = bwdist(msk);
            msk = false(size(C0));
            msk([1,end], :) = 1;
            D2 = bwdist(msk);
            D = (D1 .* D2) .^ 0.5;
        end
        D = ifftshift(D / max(D(:))) + 0.25;
        C0 = C0 ./ D;
    end
    bsz = [size(C0,1), size(C0,2)];
    if apply_mask % || dis_thresh < 2 * max(bsz)
        exp_dyx = expected_dyx - floor(expected_dyx ./ bsz) .* bsz;
        exp_dyx = round(exp_dyx) + 1;
        exp_dyx = min(max(exp_dyx, 1), bsz);
        maskc = zeros(size(C0), class(C0));
        indx = sub2ind(size(C0), exp_dyx(:,1), exp_dyx(:,2), (1:size(C0,3))');
        maskc(indx) = 1;
        scl_msk = min([dis_thresh/3, min(bsz)/100]);
        if scl_msk > 2
            maskc_ds = imresize(maskc, 1/scl_msk, 'bilinear');
            maskc_ds = imgaussfilt(maskc_ds, dis_thresh/scl_msk,'Padding', 'circular');
            maskc = imresize(maskc_ds, bsz, 'bilinear');
        else
            maskc = imgaussfilt(maskc, dis_thresh,'Padding', 'circular'); 
        end
        maskc = (maskc / max(maskc(:)) + 0.5)/1.5;
        C = maskc .* C0;
    else
        C = C0;
    end
    % C = real(ifft2(F./imgaussfilt(abs(F), 10)));
    maxval = max(max(C, [], 1), [], 2);
    mask_max = C == maxval;
    ind_max = find(mask_max(:));
    [idx1, idx2, idx3] = ind2sub(size(C), ind_max);
    yy = round(accumarray(idx3(:), idx1(:), [Nimg, 1], @mean));
    xx = round(accumarray(idx3(:), idx2(:), [Nimg, 1], @mean));
    zz = (1:Nimg)';
    ind_max = sub2ind(size(C), yy(:), xx(:), zz(:));
    mask = zeros(size(C), 'single');
    mask(ind_max) = 1;
    scl_msk = min([mask_disk/3, min(bsz)/100]);
    if scl_msk > 2
        mask_ds = imresize(mask, 1/scl_msk, 'bilinear');
        mask_ds = imgaussfilt(mask_ds, mask_disk/scl_msk,'Padding', 'circular');
        mask = imresize(mask_ds, bsz, 'bilinear');
    else
        mask = imgaussfilt(mask, mask_disk,'Padding', 'circular'); 
    end
    % mask = imgaussfilt(mask, mask_disk, 'Padding', 'circular');
    mask = mask/max(mask(:));
    mask = mask < exp(-0.5);
    % C0 = C;
    C = C0 .* mask;
    if flipcmp
        Ct = real(ifft2(IMG1.*IMG2));
        C = max(Ct,C);
    end
    maxval2 = max(max(C, [], 1), [], 2);
    conf = maxval(:)./max(maxval2(:), 1e-3);
    conf = 1 - 1./conf;
    dyx = [yy(:) - 1, xx(:) - 1];
    dyx = dyx - round((dyx - expected_dyx)./[imght, imgwd]) .* [imght, imgwd];
    if ~apply_mask
        dis = sum((dyx - expected_dyx).^2,2).^0.5;
        conf(dis > dis_thresh) = 0;
    end
    dyx = dyx - szshift;
end
