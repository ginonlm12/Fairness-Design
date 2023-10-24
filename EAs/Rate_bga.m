function obj_value = Rate_bga(x,AP,Ter,M,K,obj)

    alpha = reshape(x, K, M+1);
    
    N= 100; % Số lượng anttenna vệ tinh
    f = 20000; % Frequency in MHz
    fc = f/1e3; % Đơn vị GHz
    lambdaWave = 3e8/(fc*1e9); % Chuyển đổi đơn vị từ GHz sang MHz
    dH =lambdaWave; % The antenna spacing in the horizontal direction - first mentioned in Page 3
    MH =  sqrt(N); 
    B=100; % Băng thông hệ thống, đơn vị MHz
    tau = K; % Tín hiểu pilot
    tauc = 10000; % Số lượng tín hiệu/ subcarriers của một symbol

    hbs = 35; % Độ cao của APs, đơn vị m
    hut = 1.5; % Độ cao của người dùng, đơn vị m
    h=5; % Average building// chiều cao trung bình của các tòa nhà 
    W=20;% Average street width
    L1=  161.04 -7.1*log10(W) + 7.5*log10(h) - (24.37 -3.7*(h/hbs)^2)*log10(hbs) ...
        - (3.2*((log10(hut))^2) -4.97) - (43.42 - 3.1*log10(hbs))*3 ;% + (43.42 - 3.1*log10(hbs))
                % Giá trị của hệ số truyền tải không dây


    %aL = (1.1*log10(f)-0.7)*Hm-(1.56*log10(f)-0.8);
    %L = 46.3+33.9*log10(f)-13.82*log10(Hb)-aL;

    power_fmax = 10^(5/10); %Power per symbol in Watt
    GA = 5; %dBi;
    GSRX = 26.9; %dBi;
    NFA = 7; %9
    NFS = 1.2;%1.2
    Pda = power_fmax;%nomalized receive SNR
    Pds = power_fmax;%nomalized receive SNR
    Pua=Pda;
    Pus = Pds;
    StepVal = 1e19;
    sigma2a = StepVal*(10^((-203.975+10*log10(B*10^6)+NFA)/10))/(10^(2*GA/10));
    sigma2s = StepVal*(10^((-203.975+10*log10(B*10^6)+NFS)/10))/(10^((GA+GSRX)/10));
                % Liên quan đến phương sai trong phân phối của BNN độ nhiễu trắng
    NbrLargeFade = 10; 
    D=20; %Độ rộng của vùng xem xét, đơn vị km2
    % Satellite
    %psi = pi/6;
    RE = 6371e3; % Bán kính của trái đất
    h0 = 400e3; % Độ cao của vệ tinh

    %Shadow fading
    sigma_shda = 8; %in dB
    sigma_shds = 2.7; %in dB
    xi=((15*pi)/180)^2; %15 degree
    sigmav1 =xi;
    RateMCBothSASave = zeros(NbrLargeFade,K);
    RateCFBothSASave = zeros(NbrLargeFade,K);
    RateMCSSave = zeros(NbrLargeFade,K);
    RateCFSSave = zeros(NbrLargeFade,K);
    RateMCASave = zeros(NbrLargeFade,K);
    RateCFASave = zeros(NbrLargeFade,K);
    %AP = csvread('AP.csv');
    %Ter = csvread('Ter.csv');
    %{
    BETAA = zeros(M,K);
    dist=zeros(M,K);
    for m=1:M
            for k=1:K
            dist(m,k) = norm(AP(m,:)-Ter(k,:)); % Tính toán khoảng cách

           % if dist(m,k)<d0
            %     betadB=-L - 35*log10(d1) + 20*log10(d1) - 20*log10(d0);
            %elseif ((dist(m,k)>=d0) && (dist(m,k)<=d1))
            %     betadB= -L - 35*log10(d1) + 20*log10(d1) - 20*log10(dist(m,k));
            %else
            %    betadB = -L - 35*log10(dist(m,k)) + sigma_shda*randn(1,1); %large-scale in dB
            %end
                %betadB = -13.54 - 39.08*log10(dist(m,k)*1e3) -20*log10(fc) + 0.6*(Hm-1.5) + sigma_shda*randn(1,1);
                betadB = -L1 -   20*log10(fc) - (43.42 - 3.1*log10(hbs))*(log10(dist(m,k)*1e3));
                BETAA(m,k)= StepVal*10^(betadB/10);
                    % Tính toán hệ số large scale fading giữa AP m và user k
                    % Đổi đơn vị dB sang đơn vị tuyến tính

            end
    end
    %}
    for IterLarge =1:1
        %IterLarge
        %AP=unifrnd(-D/2,D/2,M,2);
        % disp(AP);
        % fprintf("------\n");
        %Ter=unifrnd(-D/2,D/2,K,2);
        % disp(Ter);
                % Random tọa độ APs và người dùng trong theo phân phối chuẩn
        
        BETAA = zeros(M,K);
        dist=zeros(M,K);
        % Mục đích: Tính hệ số large scale fading giữa AP m và user k
        for m=1:M
            for k=1:K
            dist(m,k) = norm(AP(m,:)-Ter(k,:)); % Tính toán khoảng cách

           % if dist(m,k)<d0
            %     betadB=-L - 35*log10(d1) + 20*log10(d1) - 20*log10(d0);
            %elseif ((dist(m,k)>=d0) && (dist(m,k)<=d1))
            %     betadB= -L - 35*log10(d1) + 20*log10(d1) - 20*log10(dist(m,k));
            %else
            %    betadB = -L - 35*log10(dist(m,k)) + sigma_shda*randn(1,1); %large-scale in dB
            %end
                %betadB = -13.54 - 39.08*log10(dist(m,k)*1e3) -20*log10(fc) + 0.6*(Hm-1.5) + sigma_shda*randn(1,1);
                betadB = -L1 -   20*log10(fc) - (43.42 - 3.1*log10(hbs))*(log10(dist(m,k)*1e3));
                BETAA(m,k)= StepVal*10^(betadB/10);
                    % Tính toán hệ số large scale fading giữa AP m và user k
                    % Đổi đơn vị dB sang đơn vị tuyến tính

            end
        end

        BETASToTal = zeros(K,1);
        BETASNonLoS = zeros(K,1);
        RAll = zeros(N,N,K);
        RAllSqrt = zeros(N,N,K);
        KRician = 6.8*ones(K,1); %Table 6.7.2-1a since psik approximate 90 degree % Hệ số của phân phối Rician
        LoSChannels = zeros(K,N);

        % Mục đích: Tính Rk - Page 3
        for k =1:K
            [Azimuk,Elek,~] = cart2sph(1e3*Ter(k,1)-300e3,1e3*Ter(k,2)-300e3,h0); %Satellite Location (300e3,300,h0)
                % Tính góc phương vị và góc độ cao giữa satellite và user k
            d = sqrt(((RE^2)*(sin(Elek))^2) + h0^2 + 2*h0*RE) - RE*sin(Elek);
                % Tính khoảng cách giữa satellite và user k - Page 11
            
            betasdB = -32.45 - 20*log10(fc) - 20*log10(d);% +sigma_shds*randn(1,1);
            % betasdB = -32.45 - 20*log10(fc) - 20*log10(d) +sigma_shds*T;
                % Tính hệ số large-scale fading giữa satellite và user k
            phik = (pi/32);%*rand(1);
            % phik = (pi/32)*TT;
            
            temp = 10*pi*sin(phik);
            Val = 4*abs(besselj(1,temp)/temp)^2;
            BETASToTal(k) = StepVal*Val*(10^(betasdB/10));
            BETASNonLoS(k)= BETASToTal(k)/(KRician(k) +1); % Page 3

            kSvec = ([cos(Azimuk)*cos(Elek), sin(Azimuk)*cos(Elek), sin(Elek)]')*2*pi/lambdaWave;
                   % Tính wave form vector l(θk, ωk) - Page 3
            hbarrdk = zeros(N,1);
            for n = 1:N
                    um = [0; mod(n-1, MH)*dH; floor((n-1)/MH)*dH];
                    hbarrdk(n) = exp(1i*kSvec'*um);
            end
            LoSChannels(k,:) =  sqrt((KRician(k)*BETASToTal(k))/(KRician(k)+1))*hbarrdk;
                    % Tính thành phân LoS - Page 3
            Raz = zeros(MH,MH);
            Rel = zeros(MH,MH);
            for nv = 1:MH
                for nh =1:MH
                    Temp1 = 1i*(2*pi*dH/lambdaWave)*(nh -nv)*cos(Elek);
                    Temp2 = -0.5*((xi*2*pi*dH/lambdaWave)^2)*((nh -nv)^2)*(sin(Elek)^2);
                    Rel(nv,nh) = exp(Temp1)*exp(Temp2);

                    D2= (2*pi*dH/lambdaWave)*(nh -nv)*sin(Elek);
                    D3= (xi*2*pi*dH/lambdaWave)*(nh -nv)*cos(Elek);
                    sigmatilde = (sin(Elek))*sigmav1;
                    D5 = (D3^2*sigmav1^2) +1;
                    Tempv1 = - ((D3^2)*(cos(Azimuk))^2)/(2*D5);
                    Tempv2 = 1i*D2*cos(Azimuk)/D5;
                    Tempv3 = -0.5*((D2*sigmatilde)^2)/D5;
                    Raz(nv,nh) = exp(Tempv1)*exp(Tempv2)*exp(Tempv3)/sqrt(D5);
                end
            end
            Rk = BETASNonLoS(k)*kron(Raz, Rel); % Tích kronecker
                %  Tính spatial correlation matrices of a planar antenna array
                % or Tính phương sai trong phân phối của gk
            RAll(:,:,k) =  Rk; 
            Rksqrt = sqrtm(Rk);
            RAllSqrt(:,:,k) = Rksqrt;

        end
    [~, GammaMatrixA, PhiMatrixS] =  Func_C_Gamma_Letter(M,N,K,Pua,Pus, tau, BETAA, RAll,sigma2a, sigma2s, alpha);   
    %Both Satellite and AP
    power_f = power_fmax*ones(1,K);
    modeS = 1; modeA =1;
    RateCF= Func_RateClosedForm_Letter_(RAll,PhiMatrixS,GammaMatrixA,BETAA,LoSChannels,power_f,Pus,B,tau,tauc,K,sigma2s, sigma2a, modeS, modeA, alpha, M);
    end

    %obj_value = 0;
    %obj_value_0alpha = 0;
    %{
    for k = 1:K
        obj_value = obj_value +  RateCF(k);
        obj_value_0alpha = obj_value_0alpha +  RateCFwithoutalpha(k);   
    end
    %}
    %obj_value = sum(RateCF(:));
    % obj_value_0alpha = sum(RateCFwithoutalpha(:)); 
    obj_value = obj(RateCF(:));
    %obj_value_0alpha = prod(RateCFwithoutalpha);   
    %obj_value = min(RateCF);
    %obj_value_0alpha = min(RateCFwithoutalpha);  
    
    
    %fprintf("alpha = 1 ----------------:");
    %disp(RateCFwithoutalpha');
    %obj_value_0alpha = min(RateCFwithoutalpha);
    %disp(obj_value_0alpha);
    %fprintf("--------------------------:");
    %disp(RateCF');
    %disp(obj_value);
    %obj_value = min(RateCF);
    
end
