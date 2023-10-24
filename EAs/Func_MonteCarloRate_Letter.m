function RateMC=Func_MonteCarloRate_Letter(GMatrixS, GMatrixA, GEstMatrixS,GEstMatrixA,ULNoiseS, ULNoiseA, power_f,B,tau,tauc,K,modeS,modeA)
    % Mục đích: Cập nhật ergoic throughput thông qua ma trận SINK_k dựa theo formula (15)&(16) 
    RateMC = zeros(K,1);
     for k =1:K 
        ukconj = conj(squeeze(GEstMatrixS(:,k,:))); 
        gk = squeeze(GMatrixS(:,k,:)); 
        umkconj = conj(squeeze(GEstMatrixA(:,k,:))); 
        gmk =  squeeze(GMatrixA(:,k,:)); 
        zkk = modeS * sum(ukconj.*gk,1) + modeA * sum(umkconj.*gmk,1);  % 13 
        
        Numk = power_f * abs(mean(zkk))^2; 
        Denomk = modeS * mean(abs(sum(ukconj.*ULNoiseS)).^2); 
        Denomk = Denomk + modeA * mean(abs(sum(umkconj.*ULNoiseA)).^2); 
        for kprime =1:K 
            gkprime = squeeze(GMatrixS(:,kprime,:));
            gmkprime =  squeeze(GMatrixA(:,kprime,:));
            zkkprime = modeS*sum(ukconj.*gkprime,1) + modeA*sum(umkconj.*gmkprime,1);
            Temp = power_f*mean(abs(zkkprime).^2);
            Denomk = Denomk + Temp;
        end
        Denomk = Denomk - Numk;
        
        RateMC(k) = B*(1-tau/tauc)*log2(1 + Numk/Denomk);
    end
end