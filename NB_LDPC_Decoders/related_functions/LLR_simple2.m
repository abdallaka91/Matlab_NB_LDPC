% I = LL(b) = 2*y/sigm^2
% P(b=0/y) = exp(I)/(exp(I)+1)
% P(b=1/y) = 1/(exp(I)+1)
%LLR(b)=log[P(b=0/y)/P(b=1/y)]=2*y/sigm^2
function LLR = LLR_simple2(y, p, sigm)
N = size(y,1);
q = 2^p;
I = 2*y/(sigm^2);
alph = 0:q-1;
alph_bin = de2bi(alph,p);
LLR = zeros(q, N);
LLR_ = zeros(q,1);
for n = 1 : N
    for i = 1 : q
        x = alph_bin(i,:);
        tt = 0;
        for j = 1 : p
            xj = x(j);
            if xj==1
                tt = tt - I(n,j); %% if LLR(c) = log[Pr(x=0)/y)/Pr(x=c)/y] => +I(j)
            end
        end
        LLR_(i) = tt;
    end
    LLR(:,n) = LLR_; %
end


