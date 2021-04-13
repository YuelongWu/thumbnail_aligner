function [des, orientation] = orient_2nd_moment(des)
    beam_num = size(des, 1);
    proj_num = size(des, 2);
    kp_num = size(des, 3);

    RY = repmat((1:beam_num)', 1, proj_num);
    desY = nansum(RY.^2 .* des.^2, 1) ./ nansum(des.^2, 1);
    [~, orientation] = max(desY, [], 2);
    orientation = orientation - 1;
    idxY = repmat(orientation, beam_num, proj_num, 1);
    F = fft(des, proj_num, 2);
    idxK = repmat(0:(proj_num - 1), beam_num, 1);
    F = F.*exp(1i*2*pi/proj_num * repmat(idxK,1,1,kp_num).*idxY);
    des = abs(ifft(F, proj_num, 2));
end