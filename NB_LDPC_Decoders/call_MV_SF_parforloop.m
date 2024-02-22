clear;
restoredefaultpath
comput_SER_BER = false;

v0v1 = [2.4 1]; %1_0.5 gave 9e-5 FER after 1e5 frames

vi = 2;
nui = 2; % nb of locations to change
max_err_cnt = 50;
max_gen = 0.5e6;
max_iter = 10;
ebn0 = 3.0:0.2:4.6; %dB
H_matrix_mat_fl_nm = '204.102.3.6.16';


pth1 = (fullfile(pwd, 'related_functions\'));
addpath(pth1);
pth2 = (fullfile(pwd, 'related_variables\'));
pth3 = (fullfile(pwd, 'related_variables\GF_arithm'));
pth4 = (fullfile(pwd, 'related_variables\alists\'));
pth5 = (fullfile(pwd, 'related_variables\alists\matrices\'));
pth6 = (fullfile(pwd, 'results\'));
pth7 = (fullfile(pwd, 'helpful functions\'));
addpath(pth7);

load([fullfile(pth4, H_matrix_mat_fl_nm) '.mat']);
H = h;

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

lgnd = ['MV_SF__' H_matrix_mat_fl_nm '_' num2str(max_iter) '_Iter_V0_'...
    num2str(v0v1(1)) '_V1_' num2str(v0v1(2))];

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
alph1 = 1:vi;
combs = alph1;
for j0 = 2:nui
    combs = combvec(combs, alph1);
end
combs = combs';

snr_cnt = length(sigma);
FERstat = zeros(snr_cnt,1);
gen_seq_cnt = zeros(snr_cnt,1);
FER = zeros(snr_cnt,1);


parforN = 100;

for i0 = 1 : snr_cnt
    FER_ = 0;
    gen_seq_cnt_ = 0;
    msg = sprintf("EbNo = %.3f dB, V0_V1= [%.3f %.3f], Error frames/Total frames = %d/%d => FER = %.8f\n",...
        ebn0(i0), v0v1(1), v0v1(2), FER(i0), gen_seq_cnt(i0), FER(i0)/gen_seq_cnt(i0));
    fprintf(msg)
    sigm =sigma(i0);
    while FER(i0) < max_err_cnt && gen_seq_cnt(i0)<max_gen
        parfor j = 1 : parforN
            gen_seq_cnt_ = gen_seq_cnt_+1;
            nse = sigm*randn(size(y_bin));
            y_bin_nse = y_bin + nse;
            LLR_2 = LLR_BPSK_GFq_2D(y_bin_nse, sigm);
            [iter, dec_seq, success_dec] = nb_ldpc_MV_SF(LLR_2, vi,v0v1, nui, max_iter, mul_mat, add_mat, div_mat,combs, h);
            rec_info_seq = dec_seq(N+1-K:N);
            nd = sum(rec_info_seq~=info_seq);
            if nd ~=0
                FER_ = FER_ +1;
            end
        end
        gen_seq_cnt(i0) = gen_seq_cnt_;
        FER(i0) = FER_;
        FERstat(i0)=FER(i0)/gen_seq_cnt(i0);

        fprintf(repmat('\b',1,length(char(msg))));
        msg = sprintf("EbNo = %.3f dB, V0_V1= [%.3f %.3f], Error frames/Total frames = %d/%d => FER = %.8f\n",...
            ebn0(i0), v0v1(1), v0v1(2), FER(i0), gen_seq_cnt(i0), FER(i0)/gen_seq_cnt(i0));
        fprintf(msg)

        



    end
end





