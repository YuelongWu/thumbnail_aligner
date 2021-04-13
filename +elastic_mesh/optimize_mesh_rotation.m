function [Ms, R] = optimize_mesh_rotation(Ms)
    tmpvar = vertcat(Ms{:});
    tmpvar = vertcat(tmpvar.TR);
    pts0 = vertcat(tmpvar.Points);
    Nang = 180;
    thetas = linspace(-89, 90, Nang);
    imgszs = nan(Nang,2);
    minpt = nan(Nang,2);
    for t = 1:Nang
        tht = thetas(t);
        ss = sind(tht);
        cs = cosd(tht);
        R = [cs, ss; -ss, cs];
        pts = pts0 * R;
        imgszs(t,:) = range(pts);
        minpt(t, :) = min(pts);
    end
    areas = imgszs(:,1) .* imgszs(:,2);
    [~, t] = min(areas + 0.5 * (imgszs(:,1)<imgszs(:,2)));
    tht = thetas(t);
    ss = sind(tht);
    cs = cosd(tht);
    offst = minpt(t,:);
    R = [cs, ss; -ss, cs];
    for k = 1:numel(Ms)
        Ms{k}.TR.Points = Ms{k}.TR.Points * R - offst;
    end
end