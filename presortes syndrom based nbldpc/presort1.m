clear
load arith_4.mat
q = 4;
dc = 4;
nm=dc-1;
di = [0, nm-1, 2];

Upc_LLR = [0 0 0 0;
          10 12 5 9;
          15 16 7 13;
          20 25 12 17];

Upc_symb = [1 2 0 3;
    0 3 1 0;
    2 0 2 1;
    3 1 3 2];

deviation_lst = syndrome_based_list_dev(dc, di);
D_sz = size(deviation_lst,1);
synd_symb = zeros(D_sz,1);
synd_LLR = synd_symb;
x_gf_vectors = zeros(size(deviation_lst));
x_LLR_vectors = zeros(size(deviation_lst));

c = zeros(D_sz,1);

temp_symb = zeros(1, dc);
temp_LLR = temp_symb;
for i = 1 : D_sz
    for j = 1 : dc
        temp_symb(1, j) = Upc_symb(deviation_lst(i, j), j);
        temp_LLR(1, j) = Upc_LLR(deviation_lst(i, j), j);
    end
    x_gf_vectors(i,:) = temp_symb;
    x_LLR_vectors(i,:) = temp_LLR;
    synd_symb(i) = sum_arr_gf_dec(temp_symb, add_mat);
    synd_LLR(i) = sum(temp_LLR);
end

[synd_LLR,sort_i]= sort(synd_LLR);
deviation_lst = deviation_lst(sort_i,:);
synd_symb = synd_symb(sort_i);
x_gf_vectors = x_gf_vectors(sort_i,:);
x_LLR_vectors= x_LLR_vectors(sort_i,:);
