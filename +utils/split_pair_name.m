function secnames = split_pair_name(pairname)
   Np = numel(pairname);
   pairname = split(pairname, '_');
   pairname = reshape(pairname, Np, numel(pairname)/Np);
   Np = size(pairname, 2);
   if Np == 2
       secname1 = pairname(:,1);
       secname2 = pairname(:,2);
   else
       secname1 = join(pairname(:,1:(Np/2)), '_');
       secname2 = join(pairname(:,(Np/2+1):end), '_');
   end
   secnames = [secname1, secname2];
end

