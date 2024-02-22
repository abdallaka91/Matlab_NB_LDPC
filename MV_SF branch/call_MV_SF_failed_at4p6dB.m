clear;
restoredefaultpath
comput_SER_BER = false;

v0v1v2 = [ 1 0.8 0.4];

vi = 2;
nui = 2; % nb of locations to change
max_err_cnt = 100;
max_iter = 10;
ebn0 = 4.6:0.2:4.6; %dB
H_matrix_mat_fl_nm = '204.102.3.6.16';
NN0 = 100;


pth1 = (fullfile(pwd, 'related_functions\'));
addpath(pth1);
pth2 = (fullfile(pwd, 'related_variables\'));
pth3 = (fullfile(pwd, 'related_variables\GF_arithm'));
pth4 = (fullfile(pwd, 'related_variables\alists\'));
pth5 = (fullfile(pwd, 'related_variables\alists\matrices\'));
pth6 = (fullfile(pwd, 'results\'));

load([fullfile(pth4, H_matrix_mat_fl_nm) '.mat']);
H = full(h);

p = ceil(log2(max(max(H))+0.1));
q = 2^p;
words = (0:q-1);

fl_nm = ['arith_' num2str(q) '.mat'];
if  exist(fullfile(pth3, fl_nm), 'file') == 2
    load(fullfile(pth3, fl_nm));
else
    add_mat = GF_arithm_matrix(q, 'add');
    mul_mat = GF_arithm_matrix(q, 'mul');
    div_mat = GF_arithm_matrix(q, 'div');
    save(fullfile(pth3, ['arith_' num2str(q) '.mat']), 'add_mat' ,'mul_mat','div_mat')
end

M = size(H, 1);
N = size(H, 2);
K = N-M;

p1 = 1;
Rate = p1*K/N; %p1 is nb of bits per channel use with the modulation, for example for bpsk it is 1

ebn0_n = 10.^(ebn0/10);
N0 = 1./(Rate*ebn0_n);
sigma = sqrt(N0/2);
snr = -10*log10(2*sigma.^2);

%%
info_seq = 0*randi([0 q-1], 1, K);
info_seq_bit=(fliplr(de2bi(info_seq,p)))';
info_seq_bit=info_seq_bit(:);
code_seq = zeros(1,N);
y_bin0 = de2bi(code_seq,p);
y_bin = (-1).^y_bin0;
H = sparse(H);
%%
str_cn_vn = cell(M,1);
dc = zeros(M,1);
for i = 1 : M
    str_cn_vn{i, 1} = find(h(i,:));
    dc(i) = length(str_cn_vn{i});
end

alph1 = 1:vi;
combs = alph1;
for j0 = 2:nui
    combs = combvec(combs, alph1);
end
combs = combs';

Tml = cell(M,1);
R = cell(M,1);
Vmn = cell(M,1);
Wmn = cell(M,1);
Dwnm = cell(M,1);
for i = 1 : M
    Tml{i} = zeros(dc(i), size(combs,1));
    R{i} = zeros(dc(i), size(combs,1));
    Vmn{i} = zeros(dc(i),size(combs,1));
    Dwnm{i} = zeros(dc(i),q);
    Wmn{i} = zeros(dc(i),q);
end

v0v1v2_cnt = size(v0v1v2, 1);
FERstat = zeros(v0v1v2_cnt,1);


saved_failed_y = cell(100);
tic

load('failed_at_4p6dB_204_102');

v0 = (0.1:0.1:2.5);
v0v1 = combvec(v0, v0)';

v0v1v2_cnt = size(v0v1,1);
gen_seq_cnt = zeros(v0v1v2_cnt,1);
FER = zeros(v0v1v2_cnt,1);

for i0 = 1 : v0v1v2_cnt
    FER_ = 0;
    FER1 = 0;
    i1 = 1;
%     v0v1v2_
        parfor  gen_seq_cn0= 1:11
            %             nse = sigma*randn(size(y_bin));
            %             y_bin_nse = y_bin + nse;
            y_bin_nse  =saved_failed_y{gen_seq_cn0};

            LLR_2 = LLR_BPSK_GFq_2D(y_bin_nse, sigma);
            [~,rec_seq_HD]=max(LLR_2, [],1);
            rec_seq_HD = rec_seq_HD-1;
            nd1 = sum(rec_seq_HD~=code_seq);

            %             [iter, dec_seq, success_dec] = nb_ldpc_MV_SF_dist(LLR_2, vi,v0v1v2_, nui, max_iter, mul_mat, add_mat, div_mat,combs, h,...
            %                 dc, str_cn_vn, Tml, R, Vmn, Dwnm, Wmn);
            [iter, dec_seq, success_dec] = nb_ldpc_MV_SF(LLR_2, vi,v0v1(i0,:), nui, max_iter, mul_mat, add_mat, div_mat, combs, h); 
            nd2 = sum(dec_seq~=code_seq);


            rec_info_seq = dec_seq(N+1-K:N);
            nd = sum(rec_info_seq~=info_seq);
            if nd2 ~=0
                FER_ = FER_+1;
                saved_failed_i = y_bin_nse;
            end
            %         FERstat(i0)=FER(i0)/gen_seq_cnt(i0);
            %         if FER_ == max_err_cnt
            %             gen_seq_cn0 = max_gen;
            %         end
        end
        toc
        FER(i0) = FER_;
        gen_seq_cnt(i0) = gen_seq_cnt(i0)+NN0;
%         clc
        fprintf("%d: V0_V1= [%.2f %.2f], Error frames/Total frames = %d/%d => FER = %.8f\n",...
            i0, v0v1(i0, 1), v0v1(i0, 2), FER(i0), gen_seq_cnt(i0), FER(i0)/gen_seq_cnt(i0))
        if FER(i0)~=FER1
%             saved_failed_y{i1} = saved_failed_i;
%             i1 = i1 + 1;
            FER1 = FER(i0);
        end
        stem(1:i0, FER(1:i0), '.')
        pause(0.01)


end




