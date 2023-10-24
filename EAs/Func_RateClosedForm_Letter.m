function RateCF = Func_RateClosedForm_Letter(RAll,PhiMatrixS,GammaMatrixA,BETAA,LoSChannels,power_f,Pus,B,tau,tauc,K,sigma2s, sigma2a, modeS, modeA)
    RateCF = zeros(K,1);
    for k =1:K
        Thetak = squeeze(RAll(:,:,k))*squeeze(PhiMatrixS(:,:,k))*squeeze(RAll(:,:,k));
        TempNumk = modeS *(norm(LoSChannels(k,:),2)^2 + Pus*tau*trace(Thetak)) +  modeA*sum(GammaMatrixA(:,k));
        % TempNumk = modeS *alpha(k,M+1)*(norm(LoSChannels(k,:),2)^2 + Pus*tau*trace(Thetak)) +  modeA*sum(alpha(k,1:M).*GammaMatrixA(:,k));
        % Numk = power_f(k)*abs(TempNumk).^2;
        Numk = power_f(k)*abs(TempNumk)^2;
        % Denomk = (modeS *alpha(k,M+1)*((sigma2s* norm(LoSChannels(k,:),2)^2) + sigma2s*Pus*tau*abs(trace(Thetak))) + modeA*sigma2a*sum(alpha(k,1:M).*GammaMatrixA(:,k)));
        Denomk = (modeS*((sigma2s* norm(LoSChannels(k,:),2)^2) + sigma2s*Pus*tau*abs(trace(Thetak))) + modeA*sigma2a*sum(GammaMatrixA(:,k)));
        for kprime = 1:K
            if  kprime ~=k
                Denomk = Denomk +  modeS*power_f(kprime)*abs(conj((LoSChannels(kprime,:)))*LoSChannels(k,:).').^2;
                % Denomk = Denomk +  modeS*alpha(k,M+1)*alpha(kprime,M+1)*power_f(kprime)*abs(conj((LoSChannels(kprime,:)))*LoSChannels(k,:).').^2;
            end
            % Denomk = Denomk +  modeS*alpha(k,M+1)*alpha(kprime,M+1)*power_f(kprime)*Pus*tau*abs(conj((LoSChannels(kprime,:)))*Thetak*LoSChannels(kprime,:).');
            Denomk = Denomk +  modeS*power_f(kprime)*Pus*tau*abs(conj((LoSChannels(kprime,:)))*Thetak*LoSChannels(kprime,:).');
            
            % Denomk = Denomk +  modeS*alpha(k,M+1)*alpha(kprime,M+1)*power_f(kprime)*abs(conj((LoSChannels(k,:)))*squeeze(RAll(:,:,kprime))*LoSChannels(k,:).');
            Denomk = Denomk +  modeS*power_f(kprime)*abs(conj((LoSChannels(k,:)))*squeeze(RAll(:,:,kprime))*LoSChannels(k,:).');
            
            % Denomk = Denomk +  modeS*alpha(k,M+1)*alpha(kprime,M+1)*power_f(kprime)*Pus*tau*abs(trace(Thetak*squeeze(RAll(:,:,kprime))));
            Denomk = Denomk +  modeS*power_f(kprime)*Pus*tau*abs(trace(Thetak*squeeze(RAll(:,:,kprime))));
            
            % Denomk = Denomk +   modeA*power_f(kprime)*sum((alpha(k,1:M)*GammaMatrixA(:,k)).*(alpha(kprime,1:M)*BETAA(:,kprime)));
            Denomk = Denomk +   modeA*power_f(kprime)*sum(GammaMatrixA(:,k).*BETAA(:,kprime));
        end
        
        RateCF(k) = B*(1-tau/tauc)*log2(1 + Numk/Denomk);
        
    end
end