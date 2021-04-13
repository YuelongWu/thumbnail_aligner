function D = composite_moving_field(mxy0, mxy1)
    D0 = elastic_mesh.movcoord2displ(mxy1);
    mxy2 = imwarp(mxy0, D0, 'FillValues', nan);
    D =  elastic_mesh.movcoord2displ(mxy2);
end