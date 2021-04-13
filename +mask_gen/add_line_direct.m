function indx = add_line_direct(imgsz, yx, dyx, return_matrix)
    yx_bd = max(1, (sign(dyx) > 0).*imgsz);
    yx_bdx = yx_bd;
    yx_bdx(:,1) = (yx_bdx(:,2) - yx(:,2)).*dyx(:,1)./dyx(:,2) + yx(:,1);
    yx_bdy = yx_bd;
    yx_bdy(:,2) = (yx_bdy(:,1) - yx(:,1)).*dyx(:,2)./dyx(:,1) + yx(:,2);
    disx = sum((yx_bdx - yx).^2, 2);
    disy = sum((yx_bdy - yx).^2, 2);
    yxt = yx_bdx;
    yxt(disx>disy,:) = yx_bdy(disx>disy,:);
    indx = mask_gen.add_lines_enpts(imgsz, yx, yxt, return_matrix);
end