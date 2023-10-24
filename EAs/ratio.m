function percentages = ratio(x,K,M)
    binaryVector = convertToBinary(x);
    alpha = reshape(binaryVector, K, M+1);
    
    percentages = zeros(3,1);
    
    percentages(2) = sum(alpha(1:K,M+1) == 0); % Connected to APs only
    for i = 1:K
        if(sum(alpha(i,1:M) == 0) == M)
            percentages(1) = percentages(1) + 1;
        end
    end % Connected to Satellite only
    percentages(3) = K - percentages(2) - percentages(1); % Connected to both APs and Satellite
    %disp(percentages);                
end