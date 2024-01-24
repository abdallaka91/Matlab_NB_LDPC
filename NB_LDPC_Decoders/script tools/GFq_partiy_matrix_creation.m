clear

q = 64;
N_eye = 26;
M_eye = 4;
eye_size = q-1;

M = M_eye*eye_size;
N = N_eye*eye_size;



p = log2(q);

H = zeros(M, N);

I = eye(eye_size,eye_size);

mat1 = zeros(M_eye, N_eye);
i1 = randperm(N_eye);
mat1(1,:) = i1;

for  i = 2 : M_eye

    s=1;
    while s~=0
        per1 = randperm(N_eye);
        s = sum(sum(per1==mat1(1:i-1,:)));
    end
    mat1(i,:) = per1;
end



mat2 = zeros(M_eye, N_eye);

per1 = randperm(q-1, N_eye);
mat2(1,:) = per1;
for  i = 2 : M_eye

    s=1;
    while s~=0
        per1 = randperm(q-1, N_eye);
        s = sum(sum(per1==mat2(1:i-1,:)));
    end
    mat2(i,:) = per1;
end

H1 = H;

for  i = 1 : M_eye
    for j = 1 : N_eye
        c1 = circshift(eye(eye_size), mat1(i, j));
        H((i-1)*eye_size+1:i*eye_size, (j-1)*eye_size+1:j*eye_size) = c1;
        H1((i-1)*eye_size+1:i*eye_size, (j-1)*eye_size+1:j*eye_size) = c1*mat2(i,j);
    end
    
end

h = H1;

save(['C:\Users\Abdallah Abdallah\Documents\Abdallah france\PhD\Matlab\nb_ldpc\alists\Generated_' ...
    num2str(M) 'x' num2str(N) '_GF' num2str(q) '.mat'],'h')







