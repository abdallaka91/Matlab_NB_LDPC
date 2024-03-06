function deviation_lst = syndrome_based_list_dev(dc, di)


D_sz = 1 + dc*di(2) + nchoosek(dc, di(3))*di(3)^2;

deviation_lst = ones(D_sz, dc);
l = 2;
for i = 1 : length(di)-1

    comb0 = 2:di(i+1)+1;
    comb00 = comb0;
    for h = 1:i-1
        comb00 = combvec(comb00,comb0);
    end
    comb00 = comb00';
    cmb = nchoosek(1:dc, i);

    for j = 1 : size(cmb,1)
        idx = ones(1, dc);
        for k = 1 : size(comb00,1)
            idx(cmb(j,:)) = comb00(k,:);
            deviation_lst(l,:) = idx;
            l=l+1;
        end
    end
end