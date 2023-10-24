function maxVal = Bruteforce(AP,Ter,M,K)
    maxVal = -Inf;

    dim = K*(M+1);
    x = zeros(1, dim);
    
    for i = 2^dim-1:-1:0
        binary = dec2bin(i, dim);
        for j = 1:dim
            x(j) = str2double(binary(j));
        end
        
        y = Rate(x,AP,Ter,M,K);

        if y > maxVal
            maxVal = y;
            disp(maxVal);
        end
    end

    disp(maxVal);
end