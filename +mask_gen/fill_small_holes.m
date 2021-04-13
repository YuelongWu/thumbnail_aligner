function mask = fill_small_holes(mask, hole_sz)
   L = unique(mask(:));
   small_holes = false(size(mask));
   for k = 1:numel(L)
       maskt = mask == L(k);
       maskt_r = bwareaopen(maskt, hole_sz, 4);
       maskt = maskt & ~(maskt_r);
       small_holes = small_holes | maskt;
   end
   [~,idx] = bwdist(~small_holes);
   mask = mask(idx);
end
