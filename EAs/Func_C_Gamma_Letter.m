function [CMatrixA, GammaMatrixA, PhiMatrixS] =  Func_C_Gamma_Letter(M,N,K,Pua,Pus, tau, BETAA, RAll,sigma2a, sigma2s, alpha)
    % Mục đích: Cập nhật GammaMatrixA, PhiMatrixS - Formula (8) 
    CMatrixA = zeros(M,K);
    GammaMatrixA = zeros(M,K);
    PhiMatrixS = zeros(N,N,K);
    IN = eye(N); % Ma trận đơn vị
    for k = 1:K
        for m = 1: M
            Numcmk = sqrt(Pua*tau)*BETAA(m,k)*alpha(k,m); % Thieu y_pm,.... % Hoi thay
            Denocmk = Pua*tau*BETAA(m,k)*alpha(k,m) + sigma2a;
            CMatrixA(m,k) = Numcmk/Denocmk; % 7
            GammaMatrixA(m,k) = Numcmk*CMatrixA(m,k);
                % Formula (8)
        end
        
        % Compute PhiMatrix for satellite
        PhiMatrixS(:,:,k) = (inv(Pus*tau*squeeze(RAll(:,:,k))*alpha(k,M+1) + sigma2s*IN)); % In Formula (9)
    end
           
end