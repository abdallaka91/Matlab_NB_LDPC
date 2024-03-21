clear;%/home/abdallah/Downloads/NB_LDPC_Decoders/call_MV_SF_parforloop.m
comput_SER_BER = false;
nm = 4;
dc1 = [0 2];
di = cell(length(dc1),1);
di{1} = [0 0];
di{2} = [0 2 1];
% di{3} = [0 6 2];


max_err_cnt = 100;
max_gen = 1e6;
max_iter = 15;
ebn0 = [3.0 :0.2:5.2]; %dB
p = 6;
q = 2^p;
K = 48;

pth1 = (fullfile(pwd, 'related_functions'));
addpath(pth1);
pth2 = (fullfile(pwd, 'related_variables'));
pth3 = (fullfile(pwd, 'related_variables/GF_arithm'));
pth4 = (fullfile(pwd, 'related_variables/alists'));
pth5 = (fullfile(pwd, 'related_variables/alists/matrices'));
pth6 = (fullfile(pwd, 'results/'));


words = (0:q-1);

H_matrix_mat_fl_nm = 'N576_K288_GF64_non_exponen_form';
load([fullfile(pth4, H_matrix_mat_fl_nm) '.mat']);
h = full(h);

N = size(h,2);
M = size(h,1);

fl_nm = ['arith_' num2str(q) '.mat'];
if  exist(fullfile(pth3, fl_nm), 'file') == 2
    load(fullfile(pth3, fl_nm));
else
    add_mat = GF_arithm_matrix(q, 'add');
    mul_mat = GF_arithm_matrix(q, 'mul');
    div_mat = GF_arithm_matrix(q, 'div');
    save(fullfile(pth3, ['arith_' num2str(q) '.mat']), 'add_mat' ,'mul_mat','div_mat')
end

dev_lsts = cell(M,1);
dev_pos = cell(M,1);
str_cn_vn = cell(M,1);
dc = zeros(M,1);

for i = 1 : M
    str_cn_vn{i, 1} = find(h(i,:));
    dc(i) = length(str_cn_vn{i});
    dc11 = dc1;
    dc11(1) = dc(i)-sum(dc1);
    [lst_deviation_lst, lst_dev_pos, dev_lsts_i, dev_pos_i]=list_dev_reliabl(dc11, di);
    dev_lsts{i} = dev_lsts_i;
    dev_pos{i} = dev_pos_i;
end
disp(length(dev_lsts_i))

% M = size(H, 1);
% N = size(H, 2);
% K = N-M;

p1 = 1;
Rate = p1*K/N; %p1 is nb of bits per channel use with the modulation, for example for bpsk it is 1

ebn0_n = 10.^(ebn0/10);
N0 = 1./(Rate*ebn0_n);
sigma = sqrt(N0/2);
snr = -10*log10(2*sigma.^2);

%%
info_seq = 0*randi([0 q-1], 1, K);
% info_seq_bit=(fliplr(de2bi(info_seq,p)))';
info_seq_bit = fliplr(dec2bin(info_seq, p) - 48)';
info_seq_bit=info_seq_bit(:);
code_seq = zeros(1, N);
% code_seq = [0     3     0     2     1     1     3     1];
info_seq = code_seq(end-K+1:end);
% y_bin0 = de2bi(code_seq,p);
y_bin0 = fliplr(dec2bin(code_seq, p) - 48);
y_bin = (-1).^y_bin0;
alph_bin =  fliplr(dec2bin(words, p) - 48);
alph_bin_mod = (-1).^alph_bin;

%%
% vi = 3;
% nui=2;
% L = vi^nui;
% alph1 = 1:vi;
% combs = alph1;
% for j0 = 2:nui
%     % combs = combvec(combs, alph1);
%     combs = CombVec(combs, alph1);
% end
% combs = combs';


str_vn_cn = cell(N,1);
dv = zeros(N,1);
for j = 1 : N
    str_vn_cn{j, 1} = (find(h(:,j)))';
    dv(j) = length(str_vn_cn{j});
end

snr_cnt = length(sigma);
FERstat = zeros(snr_cnt,1);
gen_seq_cnt = zeros(snr_cnt,1);
FER = zeros(snr_cnt,1);


parforN = 100;

iter_cnt = 0;
for i0 = 1 : snr_cnt
    FER_ = 0;
    gen_seq_cnt_ = 0;
    msg = sprintf("EbNo = %.3f dB, Error frames/Total frames = %d/%d => FER = %.8f\n",...
        ebn0(i0), FER(i0), gen_seq_cnt(i0), FER(i0)/gen_seq_cnt(i0));
    fprintf(msg)
    sigm =sigma(i0);
    while FER(i0) < max_err_cnt && gen_seq_cnt(i0)<max_gen
        tic
        parfor j = 1 : parforN
            gen_seq_cnt_ = gen_seq_cnt_+1;
            nse = sigm*randn(size(y_bin));
            y_bin_nse = y_bin + nse;
            % LLR_2 = -(round(32*0.8*LLR_simple4(y_bin_nse, p,sigm, alph_bin, alph_bin_mod)));
            % LLR_2(LLR_2<-1000)=-1000;
            LLRfact = 1;
            unreliable_sat=-inf;
            LLR_2 = LLR_simple3(y_bin_nse, p,LLRfact , unreliable_sat);
            LLR_2 = round(LLR_2*4)*32;
            [iter, dec_seq, success_dec] = presorted_MVSF_1(LLR_2', max_iter, mul_mat, add_mat, div_mat, h,str_cn_vn, dc, ...
                dev_lsts,dev_pos,nm);
            %             LLRfact = 4;
            %             unreliable_sat = -inf;
            %             LLR_2 = LLR_simple3(y_bin_nse, p,LLRfact , unreliable_sat);
            % [iter, dec_seq, success_dec] = nb_ldpc_MV_SF_opt1(LLR_2, vi,[0.5 0.25], nui,L,...
            %     max_iter, mul_mat, add_mat, div_mat,combs, h);
            iter_cnt = iter_cnt+iter;
            rec_info_seq = dec_seq(N+1-K:N);
            nd = sum(rec_info_seq~=info_seq);
            if ~success_dec && nd ~=0
                FER_ = FER_ +1;
            end
        end
        TOC=toc;
        gen_seq_cnt(i0) = gen_seq_cnt_;
        FER(i0) = FER_;
        FERstat(i0)=FER(i0)/gen_seq_cnt(i0);

        fprintf(repmat('\b',1,length(char(msg))));
        msg = sprintf("EbNo = %.3f dB, Error frames/Total frames = %d/%d => FER = %.8f\n",...
            ebn0(i0), FER(i0), gen_seq_cnt(i0), FER(i0)/gen_seq_cnt(i0));
        fprintf(msg)

    end
end
%%
ebn00 = ebn0';
FERstat0 = FER./gen_seq_cnt;

FERstat0 = [0.905000000000000   0.830000000000000   0.570000000000000   0.436666666666667   0.267500000000000   0.126250000000000   0.054210526315789   0.017241379310345   0.004807692307692   0.000919117647059   0.000155303618574   0.000031000000000   0.000008000000000   0.000001000000000];
ebn00 = ebn0(1:end-2);

FERstat1 = [0.2 0.1 0.03 0.007 0.0015 0.00025 0.000055];
ebn01 = 4:0.1:4.6;

figure(1)


semilogy(ebn00, FERstat0,'bo:', 'LineWidth',1.2)
hold on
semilogy(ebn01, FERstat1,'ro:', 'LineWidth',1.2)
xlabel('E_b/N_0 (dB)')
ylabel('FER (Log scale)')
grid on
xlim([3 5])
ylim([10e-6 1])
title('GF(32), (837, 726), L=4')
legend({'our code, L=6, itr=10, v0=0.5,v1=0.25'; 'author code, L=6, itr=10'})
