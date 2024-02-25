% I = 2*y/sigm^2
% P(b=0/y) = exp(I)/(exp(I)+1)
% P(b=1/y) = 1/(exp(I)+1)
%LLR(b)=log[P(b=0/y)/P(b=1/y)]=2*y/sigm^2
function LLR = LLR_simple1(y, p, sigm)
N = size(y,1);
q = 2^p;
I = 2*y/(sigm^2);
p0 = exp(I)./(exp(I)+1);
p1 = 1./(exp(I)+1);
alph = 0:q-1;
alph_bin = de2bi(alph,p);
Pr = zeros(q,1);
LLR = zeros(q, N);
for n = 1 : N
    for i = 1 : q
        x = alph_bin(i,:);
        pi = 1;
        for j = 1 : p
            xj = x(j);
            if xj==1
                pi = pi * p1(n, j);
            else
                pi = pi * p0(n, j);
            end
        end
        Pr(i)=pi;
    end
    LLR(:,n) = log(Pr/Pr(1)); %it should be Pr(1)./Pr
end


