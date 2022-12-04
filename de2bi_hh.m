function binary_num = de2bi_hh(decimal_num, bin_len)
% 0 <= decimal_num <= 2^bin_len - 1
binary_num = zeros(1, bin_len);    
num_tmp = decimal_num;
for ss = (bin_len-1):-1:0
    if num_tmp >= 2^ss
       binary_num(bin_len - ss) = 1;
       num_tmp = num_tmp - 2^ss;
    end
end
end