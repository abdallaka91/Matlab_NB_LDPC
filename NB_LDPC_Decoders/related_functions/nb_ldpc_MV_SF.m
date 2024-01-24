function [iter, dec_seq, success_dec] = nb_ldpc_MV_SF(LLR_2, vi,v0,v1, nui, iter_max, mul_mat, add_mat, div_mat, h)

M = size(h,1);
N = size(h,2);

q = size(LLR_2,1);
Wmn = cell(M,1);
Dwnm = cell(M,1);
Tml = cell(M,1);
R = cell(M,1);
Vmn = cell(M,1);
Qnm_3d = cell(M,1);
Lnm_3d = cell(M,1);
combinations = cell(M,1);
str_cn_vn = cell(M,1);
mn_idx_m = cell(M,1);

dc = zeros(M,1);
nu = zeros(M,1);
v = zeros(M,1);
% Qnm = zeros(1,q);

Wn = LLR_2';

for i = 1 : M
    str_cn_vn{i, 1} = find(h(i,:));
    dc(i) = length(str_cn_vn{i});
    nu(i) = nui;
    v(i) = vi;
    alph1 = 1:v;
    combs = alph1;
    for j0 = 2:nu(i)
        combs = combvec(combs, alph1);
    end
    combs = combs';
    combinations{i} = combs;
    Tml{i} = zeros(dc(i), size(combs,1));
    R{i} = zeros(dc(i), size(combs,1));
    Vmn{i} = zeros(dc(i),size(combs,1));
end

iter = 0;

for i = 1 : M
    Dwnm{i} = zeros(dc(i),q);
    Wmn{i} = zeros(dc(i),q);
end

while iter<iter_max
    iter = iter+1;
    for i = 1 : M
        idx1 = str_cn_vn{i,1};

        Qnm_3d{i} = zeros(dc(i),q);
        Lnm_3d{i} = zeros(dc(i),v(i));
        for j = 1:dc(i)
            Dwnm{i}(j, :) = Wn(idx1(j),:) - Wmn{i}(j, :);
            temp1 = Dwnm{i}(j, :);
            [a,b] = sort(temp1,'descend');
            b0 = b - 1;
            Qnm = b0;
            Qnm_3d{i}(j, :) = Qnm;
            Lnm_3d{i}(j, :) = Qnm(1:v(i));
%             Dwnm{i}(j, :) = Dwnm{i}(j, b);%---------------------------------->>>>> hereeee
        end

        df = inf(1, q);
        mn_idx = zeros(dc(i),2);
        for j0 = 1 : dc(i)
            temp = Dwnm{i}(j0, :);
            df(2:q) = abs(temp(2:q) - temp(1));
            [mn_idx(j0,1), mn_idx(j0,2)] = min(df);
        end
        [a,b] = sort(mn_idx(:,1));
        mn_idx = mn_idx(b,:);
        mn_idx_m{i} = [b(1:nu(i)) mn_idx(1:nu(i),2)]; % b is edges mn, mn_idx(1:nu(i),2) is the nearest GF


        L1 = v(i)^nu(i);
        cmbs = combinations{i};
        idx_chng = sort(mn_idx_m{i}(:,1)');
        

        for l = 1 :L1
            temp = Lnm_3d{i};
            for j0 = 1 : nu(i)
                temp(idx_chng(j0),1) = temp(idx_chng(j0),cmbs(l,j0));
            end
            Tml{i}(:,l) = temp(:,1);
        end

        Sl = zeros(L1,1);
        for l = 1 :L1
            for j0 = 1 : dc(i)
                temp = mul_mat(Tml{i}(j0,l)+1, h(i,idx1(j0))+1); %Tml{i}(l,:).*(h_gf(i,idx1));
                Sl(l) = sum_arr_gf_dec([Sl(l) temp], add_mat);%sum_arr_gf(temp);
            end
            %             Sl(l)
            cc = zeros(1,dc(1));
            for j0 = 1 : length(idx1)
                cc(j0) = div_mat(Sl(l)+1, h(i,idx1(j0))+1);
                cc(j0) = add_mat(cc(j0)+1, Tml{i}(j0,l)+1);
            end
            R{i}(:, l) = cc;
        end

        %     Vmn{i} is zeros(M,dc(i),size(combs,1));
        for l = 1 :L1
            a1 = Tml{i}(:,l);
            b1 = (squeeze(Qnm_3d{i}(:,1)));
            eq = a1==b1;
            iiv0 = find(eq);
            iiv1 = find(~eq);
            Vmn{i}(iiv0,l) = v0;
            Vmn{i}(iiv1,l) = v1;


            ii2 = R{i}(:,l)+1;
            for j1 = 1 : dc(i)
                    Wmn{i}(j1, ii2(j1)) = Wmn{i}(j1, ii2(j1)) + Vmn{i}(j1, l);
                    Wn(idx1(j1), ii2(j1)) = Wn(idx1(j1), ii2(j1)) + Vmn{i}(j1, l);
            end

        end
    end
    [~, dec_seq] = max(Wn,[],2);
    dec_seq = dec_seq' - 1;

    synd = Inf(1,M);
    
    for j1 = 1 : M
        idx1 = str_cn_vn{j1};
        tempc = zeros(dc(j1),1);
        for j0 = 1 : dc(j1)
            tempc(j0) = mul_mat(dec_seq(idx1(j0))+1, h(j1,idx1(j0))+1);
        end
        synd(j1) = sum_arr_gf_dec(tempc, add_mat);

    end

    is_not_code = sum(synd)>0;
    success_dec = ~is_not_code;

    if success_dec
        break
    end
    
end