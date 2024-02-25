% Function to generate Gray codes for GF(32)
function grayCodes = generateGrayCodesGFN(N)
    n = N; % Number of bits
    m = 2^n; % Number of elements in GF(32)
    grayCodes = zeros(m, n); % Initialize array to store Gray codes

    for i = 0:m-1
        gray = bitxor(bitshift(i, -1), i); % Compute Gray code for i
        binary = dec2bin(gray, n); % Convert Gray code to binary
        grayCodes(i+1, :) = fliplr(binary - '0'); % Store binary representation
    end
end


