%Covariance generation 
%close all
clear all
%clc
% figure(3)
hold on

rho     = 0.3;         % power splitting ratio
alpha   = 0.5;         % time fraction for EH
g2   = 0:.1:10;
%PS_dB   = -10:2:30;             % transmit SNR = Ps/N0 in dB
PS      = 2;%10.^(PS_dB./10);
naN     = (10^(-7))*1e6;    % naN = -100 dBm, BW = 1 MHZ_N
ncN     = (10^(-6))*1e6;    % naN = -90 dBm,  BW = 1 MHZ_N
naF     = (10^(-7))*1e6;
ncF     = (10^(-6))*1e6;
epsilon = 3;                % pathloss exponent
dSF     = 10;                   % S-F distance
dSN     = 3;
dNF     = dSF - dSN;
LL       = 1e2;                  % path-loss at reference distance
lSN     = LL*dSN^-3;         %1.5%    % lambda
lSF     = LL*dSF^-3;
lNF     = LL*dNF^-3;
%
eta     =0.1;              % energy conversion coefficient
RthN    = .17;                % target data rate of User N bits/s/HZ_N
RthF    = .17;               % target data rate of User N bits/s/HZ_N
%[pN,pF] = PowerAllocation(RthN,RthF);
%
pN = (2^(2*RthN) - 1)/(2^(2*RthN+2*RthF)-1);
pF = 1-pN;

C=[];
Lamc=[];
Lamd=[];
var_v = 0.01;
alph=0.5;
Z=0;%0.6;
rng(1,'twister') 





N = 4; % Number of transmit antennas
M = 8; % number of receive antennas
K = 2; % Numer of users
t = 1; % Power normalization if required

var_v = 1;
NoumberOfTrial = 1e4; 
xbins1 = 0:0.025:2.05;
xbins2 = 2.25:0.25:10;
xbins  = [xbins1 xbins2];

tol  = 1e-6;
rng(1,'twister') % for reproducibility

alph = (0.5)*(rand(1,1e5) + 1i*rand(1,1e5));   % Known transmit correlation coefficent  intentially made it big for reproducing similar results    
alph1 = (0.5)*(rand(1,1e5) + 1i*rand(1,1e5));
alph2 = (0.5)*(rand(1,1e5) + 1i*rand(1,1e5));

%Genertion of Sigma_{k,mn}, R_k, T_k matrix
        count = 1;
        count_1 = 1;
        count_2 = 1;
        Sigma_All = [];
        Sigma_All_1 = [];
        Sigma_All_2 = [];
        
        Trace_All = [];
        Trace_All_1 = [];
        Trace_All_2 = [];
        
         for k =1:K
             for m = 1:M
                 for n = 1:M
                    for indx1 = 1:1:N
                           r(indx1) = alph(count)^(indx1-1);
                           r1(indx1) = alph1(count)^(indx1-1);
                           r2(indx1) = alph2(count)^(indx1-1);
                    end
                    count =   count + 1;
                    count_1 =   count_1 + 1;
                    count_2 =   count_2 + 1;
                    if n >= m
                     
                        Sigma{k,m,n}=toeplitz(r);
                        SigmaF{k,m,n}=toeplitz(r1);
                        SigmaNF{k,m,n}=toeplitz(r2);
                    
                    else
                            
                        Sigma{k,m,n}=Sigma{k,n,m}';
                        SigmaF{k,m,n}=SigmaF{k,n,m}';
                        SigmaNF{k,m,n}=SigmaNF{k,n,m}';
                        
                    end
                    
                    Sigma_All = [Sigma_All Sigma{k,m,n}];
                    Sigma_All_1 = [Sigma_All_1 SigmaF{k,m,n}];
                    Sigma_All_2 = [Sigma_All_2 SigmaNF{k,m,n}];
                    
                    Trace_All = [Trace_All trace(Sigma{k,m,n})];
                    Trace_All_1 = [Trace_All_1 trace(SigmaF{k,m,n})];
                    Trace_All_2 = [Trace_All_2 trace(SigmaNF{k,m,n})];
                    
                    clear r  indx1
                    clear r1  indx1
                    clear r2  indx1
                    
                 end
             end
             RT = zeros(N,N);
             RTF = zeros(N,N);
             RTNF = zeros(N,N);
             
             for indx2 = 1:M
                 RT_1 = Sigma{k,indx2,indx2};
                 RT  = RT + RT_1; 
                 
                 RT_2 = SigmaF{k,indx2,indx2};
                 RTF  = RTF + RT_2;
                 
                 RT_3 = SigmaNF{k,indx2,indx2};
                 RTNF  = RTNF + RT_3;
             end
             R_k{k} = RT;
              R_F{k} = RTF;
               R_NF{k} = RTNF;
               
             T_k{k} = K*ones(M,M);
             T_F{k} = K*ones(M,M);
             T_NF{k} = K*ones(M,M);
             
         end         
          
%Genertion of Sigma_k
    coldist1     = size(Sigma_All, 2) / K * ones(1, K); %split the number of columns by K and replicate K times.
    splitmats1   = mat2cell(Sigma_All, N, coldist1);

    for k=1:K
        ind_mat = splitmats1{k};
        coldist2     = size(ind_mat, 2) / M * ones(1, M); %split the number of columns by M and replicate M times.
        splitmats2   = mat2cell(ind_mat, N, coldist2);
        sigma_dummy = [];
        for m = 1:M
            sigma_dummy1 = splitmats2{m};
            sigma_dummy = [sigma_dummy; sigma_dummy1];        
        end
        Sigma_k{k} = sigma_dummy;
    end
%     test1= ishermitian(Sigma_k{1})
%     test2= ishermitian(Sigma_k{2})
% keybord


%Genertion of Sigma_F
    coldist1_1     = size(Sigma_All_1, 2) / K * ones(1, K); %split the number of columns by K and replicate K times.
    splitmats1_1   = mat2cell(Sigma_All_1, N, coldist1_1);

    for k=1:K
        ind_mat_1 = splitmats1_1{k};
        coldist2_1     = size(ind_mat_1, 2) / M * ones(1, M); %split the number of columns by M and replicate M times.
        splitmats2_1   = mat2cell(ind_mat_1, N, coldist2_1);
        sigma_dummy_1 = [];
        for m = 1:M
            sigma_dummy1_1 = splitmats2_1{m};
            sigma_dummy_1 = [sigma_dummy_1; sigma_dummy1_1];        
        end
        Sigma_F{k} = sigma_dummy_1;
    end
%     test1= ishermitian(Sigma_k{1})
%     test2= ishermitian(Sigma_k{2})
% keybord



%Genertion of Sigma_NF
    coldist1_2     = size(Sigma_All_2, 2) / K * ones(1, K); %split the number of columns by K and replicate K times.
    splitmats1_2   = mat2cell(Sigma_All_2, N, coldist1_2);

    for k=1:K
        ind_mat_2 = splitmats1_2{k};
        coldist2_2     = size(ind_mat_2, 2) / M * ones(1, M); %split the number of columns by M and replicate M times.
        splitmats2_2   = mat2cell(ind_mat_2, N, coldist2_2);
        sigma_dummy_2 = [];
        for m = 1:M
            sigma_dummy1_2 = splitmats2_2{m};
            sigma_dummy_2 = [sigma_dummy_2; sigma_dummy1_2];        
        end
        Sigma_NF{k} = sigma_dummy_2;
    end
%     test1= ishermitian(Sigma_k{1})
%     test2= ishermitian(Sigma_k{2})
% keybord




%%%%%%%%%%%%%%Recieve beamvectors
            for k=1:K
                u_k=sqrt(.5)*(rand(N,1)+ 1i* rand(N,1));
                vk(:,k)=(u_k/norm(u_k));
            end
            
            for k=1:K
                u_F=sqrt(.5)*(rand(N,1)+ 1i* rand(N,1));
                vF(:,k)=(u_F/norm(u_F));
            end
            
            for k=1:K
                u_NF=sqrt(.5)*(rand(N,1)+ 1i* rand(N,1));
                vNF(:,k)=(u_NF/norm(u_NF));
            end
%Genertion of Phi_k exp (4)            
    for k=1:K
       Phi_k{k} =  (kron(eye(M),vk(:,k)')*Sigma_k{k}*(kron(eye(M),vk(:,k)))).' ;
    end

    for k=1:K
       Phi_F{k} =  (kron(eye(M),vF(:,k)')*Sigma_F{k}*(kron(eye(M),vF(:,k)))).' ;
    end
    
    for k=1:K
       Phi_NF{k} =  (kron(eye(M),vNF(:,k)')*Sigma_NF{k}*(kron(eye(M),vNF(:,k)))).' ;
    end

%Genertion of Sigma_{k,mn} exp(15)
users = [1:K];
for k = 1:K
    MATA = zeros(N,N);
    for indx3 = 1:M
        for indx4 = 1:M
           MATA1 = Sigma{k,indx3,indx4};
           MATA = MATA + MATA1; 
        end
    end
    Sum_Sigma_kmn{k} = MATA;
end

for k = 1:K
    MATA_1 = zeros(N,N);
    for indx3 = 1:M
        for indx4 = 1:M
           MATA1_1 = SigmaF{k,indx3,indx4};
           MATA_1 = MATA_1 + MATA1_1; 
        end
    end
    Sum_Sigma_kmnF{k} = MATA_1;
end



for k = 1:K
    MATA_2 = zeros(N,N);
    for indx3 = 1:M
        for indx4 = 1:M
           MATA1_2 = SigmaNF{k,indx3,indx4};
           MATA_2 = MATA_2 + MATA1_2; 
        end
    end
    Sum_Sigma_kmnNF{k} = MATA_2;
end



count = 0;
count1 = 0;
count2 = 0;

while 1
for k = 1:K
    
        X2 = vk;
        X3 = vF;
        X4 = vNF;
        
        vk_k = X2(:,k);
        vF_k = X3(:,k);
        vNF_k = X4(:,k);
        
        vi = X2;
        viF = X3;
        viNF = X4;
        
        vi(:, k) = [];
        viF(:, k) = [];
        viNF(:, k) = [];
        
            usersi = users;
            usersi(k) = [];
        
            RQ_i=0;
            RQ_F=0;
            RQ_NF=0;
            for intf = 1:size(vi,2) 
    
                vi2 = vi(:,intf)'*(Sum_Sigma_kmn{usersi(intf)}*inv(R_k{usersi(intf)}))*vi(:,intf);       %sqrt(1/(M*var_v))*(kron(eye(M), (RTk*wi(:,intf)).')'*vk_B')' ;
                vi2F = viF(:,intf)'*(Sum_Sigma_kmnF{usersi(intf)}*inv(R_F{usersi(intf)}))*viF(:,intf);       %sqrt(1/(M*var_v))*(kron(eye(M), (RTk*wi(:,intf)).')'*vk_B')' ;
                vi2NF = viNF(:,intf)'*(Sum_Sigma_kmnNF{usersi(intf)}*inv(R_NF{usersi(intf)}))*viNF(:,intf);       %sqrt(1/(M*var_v))*(kron(eye(M), (RTk*wi(:,intf)).')'*vk_B')' ;

                
                %     keyboard          
                 RQ_i = RQ_i+  vi2;
                 RQ_F = RQ_F+  vi2F;
                 RQ_NF = RQ_NF+vi2NF;
            end
        Objecfun(k) =vk_k'*(inv(R_k{k})*conj(RQ_i) * Sum_Sigma_kmn{k})*vk_k
        
        [U11 V11] =eig((inv(R_k{k})*conj(RQ_i) * Sum_Sigma_kmn{k}));
        [vala,indxa] = max(real(diag(V11)));
        dom_ka = U11(:,indxa);
        vk_new(:,k) = dom_ka;
        
        
         Objecfun1(k) =vF_k'*(inv(R_F{k})*conj(RQ_F) * Sum_Sigma_kmnF{k})*vF_k
        
        [U12 V12] =eig((inv(R_F{k})*conj(RQ_F) * Sum_Sigma_kmnF{k}));
        [valaF,indxaF] = max(real(diag(V12)));
        dom_kaF = U12(:,indxaF);
        vF_new(:,k) = dom_kaF;
        
        
        
         Objecfun2(k) =vNF_k'*(inv(R_NF{k})*conj(RQ_NF) * Sum_Sigma_kmnNF{k})*vNF_k
        
        [U13 V13] =eig((inv(R_NF{k})*conj(RQ_NF) * Sum_Sigma_kmnNF{k}));
        [valaNF,indxaNF] = max(real(diag(V13)));
        dom_kaNF = U13(:,indxaNF);
        vNF_new(:,k) = dom_kaNF;
        
end


vk = vk_new;
vF = vF_new;
vNF = vNF_new;

if abs(real(max(Objecfun)) - real(min(Objecfun))) >=tol
    abs(norm(max(Objecfun)) - norm(min(Objecfun)))
    count = count + 1;
    if count ==10
        break
    end
    else
        break
end

if abs(real(max(Objecfun1)) - real(min(Objecfun1))) >=tol
    abs(norm(max(Objecfun1)) - norm(min(Objecfun1)))
    count1 = count1 + 1;
    if count1 ==10
        break
    end
    else
        break
end


if abs(real(max(Objecfun2)) - real(min(Objecfun2))) >=tol
    abs(norm(max(Objecfun2)) - norm(min(Objecfun2)))
    count2 = count2 + 1;
    if count2 ==10
        break
    end
    else
        break
end


end

%Genertion of exp(6) again using new vk
          
    for k=1:K
       Phi_k{k} =  (kron(eye(M),vk(:,k)')*Sigma_k{k}*(kron(eye(M),vk(:,k)))).' ;
       Phi_F{k} =  (kron(eye(M),vF(:,k)')*Sigma_F{k}*(kron(eye(M),vF(:,k)))).' ;
       Phi_NF{k} =  (kron(eye(M),vNF(:,k)')*Sigma_NF{k}*(kron(eye(M),vNF(:,k)))).' ;
    end

    
    
% % % % % % %     Transmit beamvectors Random R-BF
%         wk=[];         
%         for k=1:K
%             u=sqrt(.5)*(rand(M,1)+ 1i*rand(M,1));
%             wk(:,k)=(u/norm(u));
%         end
    
    
% % % % % % %     Transmit beamvectors (PEV) 
    wk_RQ=[];
    for k = 1:K
        [U1 V1] = eig(Phi_k{k});
        [valmf,indxmf] = max(real(diag(V1)));
        dom_kmf = U1(:,indxmf);
        wk_RQ = [wk_RQ dom_kmf];
    end

% % % % % % % %     Transmit beamvectors (PEV) G-ZF (2-users only)

    wk_ZF=[];
    for k = 1:K
        [U1zf{k} V1zf{k}] = eig(Phi_k{k});
    end
    U_k1 = U1zf{1};
    U_k2 = U1zf{2};
    [eigveck1 eigvalk1] = eig(U_k2(:,end)'*Phi_k{1}*U_k2(:,end));
    [val3,indx3] = max(real(diag(eigvalk1)));
    dom_k3 = eigveck1(:,indx3);
    wk1 = U_k2(:,end)*dom_k3;
    
    [eigveck2 eigvalk2] = eig(U_k1(:,end)'*Phi_k{2}*U_k1(:,end));
    [val4,indx4] = max(real(diag(eigvalk2)));
    dom_k4 = eigveck2(:,indx4);
    wk2 = U_k1(:,end)*dom_k4; 
    wk_ZF = [wk1 wk2];

% % % % % %     %     Transmit beamvectors  (2-users only)
    %%%For 2 users
    wk_MF=[];
    [eigveck1 eigvalk1] = eig(inv(Phi_k{2} + var_v*eye(M))*Phi_k{1});
    [val1,indx1] = max(real(diag((eigvalk1))));
    dom_k1 = eigveck1(:,indx1);
    wk1 = dom_k1;
    
    [eigveck2 eigvalk2] = eig(inv(Phi_k{1} + var_v*eye(M))*Phi_k{2});
    [val2,indx2] = max(real(diag((eigvalk2))));
    dom_k2 = eigveck2(:,indx2);
    wk2 = dom_k2;
    
    wk_MF = [wk1 wk2];
%     
% %     %%%For 4 users
%     wk_RQ=[];
%     
%     [eigveck1 eigvalk1] = eig(inv(Phi_k{2}+Phi_k{3}+Phi_k{4} + var_v*eye(M))*Phi_k{1});     
%     [val1,indx1] = max(real(diag((eigvalk1))));
%     dom_k1 = eigveck1(:,indx1);
%     wk1 = dom_k1;
%     
%     [eigveck2 eigvalk2] = eig(inv(Phi_k{1}+Phi_k{3}+Phi_k{4} + var_v*eye(M))*Phi_k{2});
%     [val2,indx2] = max(real(diag((eigvalk2))));
%     dom_k2 = eigveck2(:,indx2);
%     wk2 = dom_k2;
%     
%     [eigveck3 eigvalk3] = eig(inv(Phi_k{2}+Phi_k{1}+Phi_k{4} + var_v*eye(M))*Phi_k{3});
%     [val3,indx3] = max(real(diag((eigvalk3))));
%     dom_k3 = eigveck3(:,indx3);
%     wk3 = dom_k3;
%     
%     [eigveck4 eigvalk4] = eig(inv(Phi_k{2}+Phi_k{3}+Phi_k{1} + var_v*eye(M))*Phi_k{4});
%     [val4,indx4] = max(real(diag((eigvalk4))));
%     dom_k4 = eigveck4(:,indx4);
%     wk4 = dom_k4;
%     
%     wk_RQ = [wk1 wk2 wk1 wk2];
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%      G-MF  %%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%% Simulation Setup

[A, B, Y, CDF] = Coovariance_Model_Sim_4(wk_MF, M, N, K, var_v, Phi_k,Phi_F,Phi_NF,PS,g2,alpha,rho,RthF,naN ,ncN,pN,pF,ncF,naF,eta,NoumberOfTrial,xbins);

CDF_SIM_MF = CDF;

%CDF_SIM_MF = 1/K*sum(CDF.')';

%%%%%%%%%%%%%% Analytical Part

  CDFA_MF =   Coovariance_Model_Analy_4(wk_MF, M, N, K, var_v, Phi_k,PS,g2,alpha,rho,RthF,pF,pN, tol);
    
   CDF_Anal_MF = CDFA_MF;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%      G-ZF  %%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%Analytical Part

    CDFA_ZF =   Coovariance_Model_Analy_4(wk_ZF, M, N, K, var_v, Phi_k,PS,g2,alpha,rho,RthF,pF,pN, tol);
    
   CDF_Anal_ZF = CDFA_ZF;

 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%      G-RQ(SLNR)  %%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Simulation Setup

[A, B, Y, CDF_RQ] = Coovariance_Model_Sim_4(wk_RQ, M, N, K, var_v, Phi_k,Phi_F,Phi_NF,PS,g2,alpha,rho,RthF,naN ,ncN,pN,pF,ncF,naF,eta,NoumberOfTrial,xbins);

CDF_SIM_RQ = CDF_RQ;

%%%%%%%%%%%%%%%%Analytical Part

   CDFA_RQ =   Coovariance_Model_Analy_4(wk_RQ, M, N, K, var_v, Phi_k,PS,g2,alpha,rho,RthF,pF,pN, tol);
    
   CDF_Anal_RQ = CDFA_RQ;
    
    
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%%%%%%%%%%%%%%%%%%%      Optimization Part  %%%%%%%%%%%%%%%%%%%%%%
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% jj = length(y);
% 
% %%%%%%%%%%%%%%%%%%%% Transmit Adaptive BF Only %%%%%%%%%%%%%%%%%%%%
% 
% CDFOpt1 = [];
% tic
% 
% for gammaop = 1:jj
% 
%     Transmit_BF_Only = sprintf('Process %d of  %d.',gammaop,jj)
%     ynew = y(gammaop); 
% %     Yorg= CDFA(gammaop,:);
%     
%     
%         s1 = CDF_Anal_MF(gammaop); s2 = CDF_Anal_ZF(gammaop); s3 = CDF_Anal_RQ(gammaop);
%         [val, indx] = min([s1 s2 s3]); %s2 s3
%     
%         if indx == 1
%             CDFA_fb = CDFA_MF(gammaop,:);
% 
%         elseif indx == 2
%             CDFA_fb = CDFA_ZF(gammaop,:);
%         elseif indx == 3
%             CDFA_fb = CDFA_RQ(gammaop,:);
%         end
%     CDFOpt1 = [CDFOpt1; CDFA_fb];
% end
% Ob1Time = toc
% 
% CDFSLA    = (1/K)*CDFOpt1;
% ObjFunc_thres = sum(CDFSLA')';
% 
% % %%%%%Simulation Plot for the first user

% plot(y,CDF_SIM_MF,'or');

% figure(3)
hold on
xlabel('\gamma')
ylabel('Pr( SINR < \gamma )')
plot(g2,CDF_Anal_ZF,'-o');
%plot(y(1:1:end),CDF_Anal_MF(1:1:end),'-');  
%plot(g2,CDF_SIM_MF,'-o');
plot(g2,CDF_Anal_MF,'-s');
%plot(y(1:1:end),CDF_SIM_MF(1:1:end),'*'); PS_dB 
% plot(1:5:numel(y), CDF_Anal_MF(1:5:end), 'bo')

%plot(g2,CDF_Anal_ZF,'*');

%plot(y(1:1:end),CDF_Anal_ZF(1:1:end),'-');
%plot(y(1:1:end),CDF_Anal_ZF(1:1:end),'o');


%plot(g2,CDF_SIM_RQ,'-<');
plot(g2,CDF_Anal_RQ,'->');


ylim([0 1])


