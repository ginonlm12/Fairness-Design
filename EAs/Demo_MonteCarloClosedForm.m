% Demo Closed-form Monte-Carlo
clc;clear; close all;
rng('default');
M = 40; %Number of APs
K = 20; % Number of users
N= 100; %Number of Satellite antennas

%%
% alpha = randi([0, 1], K, M+1);
% alpha = zeros(K,M+1);
  alpha = ones(K,M+1);
%%

f = 20000; % Frequency in MHz
fc = f/1e3; %in GHz
lambdaWave = 3e8/(fc*1e9);
dH =lambdaWave;
MH =  sqrt(N);
B=100; %Mhz
tau =K; %Number of pilots;
tauc = 10000;

hbs = 35; % AP height in m
hut = 1.5; % user height in m
h=5; %Average building
W=20;%Average street width
L1=  161.04 -7.1*log10(W) + 7.5*log10(h) - (24.37 -3.7*(h/hbs)^2)*log10(hbs) ...
    - (3.2*((log10(hut))^2) -4.97) - (43.42 - 3.1*log10(hbs))*3 ;% + (43.42 - 3.1*log10(hbs))

%aL = (1.1*log10(f)-0.7)*Hm-(1.56*log10(f)-0.8);
%L = 46.3+33.9*log10(f)-13.82*log10(Hb)-aL;

power_fmax = 10^(5/10); %Power per symbol in Watt
GA = 5; %dBi;
GSRX = 26.9; %dBi;
NFA = 7; %9
NFS = 1.2;%1.2
Pda = power_fmax; % nomalized receive SNR
Pds = power_fmax; % nomalized receive SNR
Pua=Pda;
Pus = Pds;
StepVal = 1e19;
sigma2a = StepVal*(10^((-203.975+10*log10(B*10^6)+NFA)/10))/(10^(2*GA/10));
sigma2s = StepVal*(10^((-203.975+10*log10(B*10^6)+NFS)/10))/(10^((GA+GSRX)/10));
NbrLargeFade = 500;
NbrSmallFade = 2e3;
D=20; %km2
% Satellite
%psi = pi/6;
RE = 6371e3;
h0 = 400e3; %Satellite's attitude

%Shadow fading
sigma_shda=8; %in dB
sigma_shds=2.7; %in dB
xi=((15*pi)/180)^2; %15 degree
sigmav1 =xi;
GMatrixS = zeros(N,K,NbrSmallFade);
RateMCBothSASave = zeros(NbrLargeFade,K);
RateCFBothSASave = zeros(NbrLargeFade,K);
RateMCSSave = zeros(NbrLargeFade,K);
RateCFSSave = zeros(NbrLargeFade,K);
RateMCASave = zeros(NbrLargeFade,K);
RateCFASave = zeros(NbrLargeFade,K);
for IterLarge =1:NbrLargeFade
    IterLarge
    AP=unifrnd(-D/2,D/2,M,2);
    Ter=unifrnd(-D/2,D/2,K,2);
    BETAA = zeros(M,K);
    dist=zeros(M,K);
    for m=1:M
        for k=1:K
        dist(m,k) = norm(AP(m,:)-Ter(k,:));

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
        end
    end
    BETASToTal = zeros(K,1);
    BETASNonLoS = zeros(K,1);
    RAll = zeros(N,N,K);
    RAllSqrt = zeros(N,N,K);
    KRician = 6.8*ones(K,1);%Table 6.7.2-1a since psik approximate 90 degree
    LoSChannels = zeros(K,N);
    for k =1:K
        [Azimuk,Elek,~] = cart2sph(1e3*Ter(k,1)-300e3,1e3*Ter(k,2)-300e3,h0); %Satellite Location (300e3,300,h0)
        d = sqrt(((RE^2)*(sin(Elek))^2) + h0^2 + 2*h0*RE) - RE*sin(Elek);
        betasdB = -32.45 - 20*log10(fc) - 20*log10(d) +sigma_shds*randn(1,1);
        phik = (pi/32)*rand(1);
        temp = 10*pi*sin(phik);
        Val = 4*abs(besselj(1,temp)/temp)^2;
        BETASToTal(k) = StepVal*Val*(10^(betasdB/10)); 
        BETASNonLoS(k)= BETASToTal(k)/(KRician(k) +1); % 4

        kSvec = ([cos(Azimuk)*cos(Elek), sin(Azimuk)*cos(Elek), sin(Elek)]')*2*pi/lambdaWave; %2
        hbarrdk = zeros(N,1);
        for n = 1:N
                um = [0; mod(n-1, MH)*dH; floor((n-1)/MH)*dH];
                hbarrdk(n) = exp(1i*kSvec'*um); %Trong công thức số 1
        end
        LoSChannels(k,:) =  sqrt((KRician(k)*BETASToTal(k))/(KRician(k)+1))*hbarrdk; %1
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
        Rk = BETASNonLoS(k)*kron(Raz, Rel); %4
        RAll(:,:,k) =  Rk; 
        Rksqrt = sqrtm(Rk); % Căn bậc 2 của 4
        RAllSqrt(:,:,k) = Rksqrt;
            % Random the channel gk between the UL transmitter of user k and the satellite receiver
        for iterloop = 1:NbrSmallFade
            gk = (LoSChannels(k,:).') +  Rksqrt*((randn(N,1) + 1i*randn(N,1))/sqrt(2));
            GMatrixS(:,k,iterloop)= gk; %The channel between the UL transmitter of user k and the satellite receiver
        end
    end
    
    GMatrixA = repmat(BETAA.^0.5,[1,1,NbrSmallFade]).*((randn(M,K,NbrSmallFade) + 1i*randn(M,K,NbrSmallFade))/sqrt(2));
    [CMatrixA, GammaMatrixA, PhiMatrixS] =  Func_C_Gamma_Letter(M,N,K,Pua,Pus, tau, BETAA, RAll,sigma2a, sigma2s);
    % Random độ nhiễu trắng tuân theo phân phối CN(0,sigma2a)
    TrainNoiseA = sqrt(sigma2a)*(randn(M,K,NbrSmallFade) + 1i*randn(M,K,NbrSmallFade))/sqrt(2); %Training Noise at APs
    TrainNoiseS = sqrt(sigma2s)*(randn(N,K,NbrSmallFade) + 1i*randn(N,K,NbrSmallFade))/sqrt(2); %Training Noise at Satellite
    
    % Compute received training signals at APs
    YpMatrixA = zeros(M,K,NbrSmallFade);
    GEstMatrixA = zeros(M,K,NbrSmallFade);
               
    for m = 1:M
        for k =1:K
            ypmk = sqrt(Pua*tau)*squeeze(GMatrixA(m,k,:)) + squeeze(TrainNoiseA(m,k,:));
                % Formula 5
            YpMatrixA(m,k,:) = ypmk;
            GEstMatrixA(m,k,:) = CMatrixA(m,k)*ypmk;
        end
    end
    
    % Compute received training signals at satellite
    YpMatrixS = zeros(N,K,NbrSmallFade);
    GEstMatrixS = zeros(N,K,NbrSmallFade);
    for k =1:K
        % Formula (6)
        ypmk = (sqrt(Pus*tau)*squeeze(GMatrixS(:,k,:)) + squeeze(TrainNoiseS(:,k,:)));
        %YpMatrixS(:,k,:) = ypmk;
        ypmkbar = sqrt(Pus*tau)*LoSChannels(k,:).';
        GEstMatrixS(:,k,:) = repmat(LoSChannels(k,:).',[1, NbrSmallFade]) + sqrt(Pus*tau)*squeeze(RAll(:,:,k))*squeeze(PhiMatrixS(:,:,k))*(ypmk - repmat(ypmkbar, [1, NbrSmallFade]));
    end
            % Formula (11)
    ULNoiseA = sqrt(sigma2a)*(randn(M,NbrSmallFade) + 1i*randn(M,NbrSmallFade))/sqrt(2); %Uplink noise AP
    ULNoiseS= sqrt(sigma2s)*(randn(N,NbrSmallFade) + 1i*randn(N,NbrSmallFade))/sqrt(2); %Uplink noise Satellite
   
   % Both Satellite and AP
   power_f = power_fmax*ones(1,K);
   modeS = 1; modeA =1;
   RateMCBothSASave(IterLarge,:)=Func_MonteCarloRate_Letter(GMatrixS, GMatrixA, GEstMatrixS,GEstMatrixA,ULNoiseS, ULNoiseA, power_f,B,tau,tauc,K,modeS,modeA);
   RateCFBothSASave(IterLarge,:) = Func_RateClosedForm_Letter(RAll,PhiMatrixS,GammaMatrixA,BETAA,LoSChannels,power_f,Pus,B,tau,tauc,K,sigma2s, sigma2a, modeS, modeA);
    
    % Satellite Only
    modeS = 1; modeA =0;
    RateMCSSave(IterLarge,:)=Func_MonteCarloRate_Letter(GMatrixS, GMatrixA, GEstMatrixS,GEstMatrixA,ULNoiseS, ULNoiseA, power_f,B,tau,tauc,K,modeS,modeA);
    RateCFSSave(IterLarge,:) = Func_RateClosedForm_Letter(RAll,PhiMatrixS,GammaMatrixA,BETAA,LoSChannels,power_f,Pus,B,tau,tauc,K,sigma2s, sigma2a, modeS, modeA);
    
    % APs Only
    modeS = 0; modeA =1;
    RateMCASave(IterLarge,:)=Func_MonteCarloRate_Letter(GMatrixS, GMatrixA, GEstMatrixS,GEstMatrixA,ULNoiseS, ULNoiseA, power_f,B,tau,tauc,K,modeS,modeA);
    RateCFASave(IterLarge,:) = Func_RateClosedForm_Letter(RAll,PhiMatrixS,GammaMatrixA,BETAA,LoSChannels,power_f,Pus,B,tau,tauc,K,sigma2s, sigma2a, modeS, modeA);
    a=5;
end
StepSize = 100;
yaxis = linspace(0,1,NbrLargeFade);
figure;
plot(sort(sum(RateMCBothSASave,2),'ascend'),yaxis,'r-', 'LineWidth',1);
hold on
plot(sort(sum(RateCFBothSASave,2),'ascend'),yaxis,'ro','MarkerIndices', 1:StepSize:NbrLargeFade);
hold on
plot(sort(sum(RateMCSSave,2),'ascend'),yaxis,'k-.', 'LineWidth',1);
hold on
plot(sort(sum(RateCFSSave,2),'ascend'),yaxis,'k+','MarkerIndices', 1:StepSize:NbrLargeFade);
plot(sort(sum(RateMCASave,2),'ascend'),yaxis,'b--', 'LineWidth',1);
hold on
plot(sort(sum(RateCFASave,2),'ascend'),yaxis,'b*','MarkerIndices', 1:StepSize:NbrLargeFade);
legend('Monte Carlo, Sat.-Terrest.','Closed-Form, Sat.-Terrest.','Monte Carlo, Sat.','Closed-Form, Sat.','Monte Carlo, Terrest.','Closed-Form, Terrest.');
xlabel('Sum Data Throughput [Mbps]');
ylabel('CDF');
grid on; box on;

figure;
plot(sort(min(RateMCBothSASave,[],2),'ascend'),yaxis,'r-', 'LineWidth',1);
hold on
plot(sort(min(RateCFBothSASave,[],2),'ascend'),yaxis,'ro','MarkerIndices', 1:StepSize:NbrLargeFade);
plot(sort(min(RateMCSSave,[],2),'ascend'),yaxis,'k-.', 'LineWidth',1);
hold on
plot(sort(min(RateCFSSave,[],2),'ascend'),yaxis,'k+','MarkerIndices', 1:StepSize:NbrLargeFade);
plot(sort(min(RateMCASave,[],2),'ascend'),yaxis,'b--', 'LineWidth',1);
hold on
plot(sort(min(RateCFASave,[],2),'ascend'),yaxis,'b*','MarkerIndices', 1:StepSize:NbrLargeFade);
legend('Monte Carlo, Sat.-Terrest.','Closed-Form, Sat.-Terrest.','Monte Carlo, Sat.','Closed-Form, Sat.','Monte Carlo, Terrest.','Closed-Form, Terrest.');
xlabel('Minimum Ergodic Data Throughput [Mbps]');
ylabel('CDF');
grid on; box on;





