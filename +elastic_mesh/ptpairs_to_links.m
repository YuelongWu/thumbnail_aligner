function links = ptpairs_to_links(ptpairs, M1, M2)
    [B1, ID1, oob1] = elastic_mesh.cart2bary(M1.TR0, ptpairs.yx1);
    rid1 = double(ptpairs.region_id(:,1));
    ID1(rid1 > 0) = ID1(rid1 > 0) + (rid1(rid1 > 0) - 1) * M1.tri_num0;
    [B2, ID2, oob2] = elastic_mesh.cart2bary(M2.TR0, ptpairs.yx2);
    rid2 = double(ptpairs.region_id(:,2));
    ID2(rid2 > 0) = ID2(rid2 > 0) + (rid2(rid2 > 0) - 1) * M2.tri_num0;
    h1 = 1 - 0.5 * oob1(:);
    h2 = 1 - 0.5 * oob2(:);
    conf = ptpairs.conf;
    conf = conf/max(1, max(conf(:)));
    mass1 = M1.mass(rid1);
    mass2 = M2.mass(rid2);
    rmass1 = 2 * mass2 ./ (mass1 + mass2);
    rmass2 = 2 * mass1 ./ (mass1 + mass2);
    links = table(B1, ID1, h1, rmass1, B2, ID2, h2, rmass2, conf);
end