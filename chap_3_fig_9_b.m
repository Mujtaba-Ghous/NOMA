close all
clear all
%
%% Simulation parameters
%
N = 3;   
K=2;% # of antenna
M=1;
rho     = .1:.1:.9;         % power splitting ratio
alpha   = .1:.1:.9;         % time fraction for EH
PS_dB   = 0;                % transmit SNR = Ps/N0 in dB
PS      = 10.^(PS_dB./10);
naN     = (10^(-7))*1e6;    % naN = -100 dBm, BW = 1 MHz
ncN     = (10^(-6))*1e6;    % naN = -90 dBm,  BW = 1 MHz
naF     = (10^(-7))*1e6;
ncF     = (10^(-6))*1e6;
epsilon = 3;                % pathloss exponent
dSF = 10;                   % S-F distance
dSN = 3;
dNF = dSF - dSN;
L   = 1e3;                  % path-loss at reference distance
%
lSN = L*dSN^-3;             % lambda
lSF = L*dSF^-3;
lNF = L*dNF^-3;
%
eta     = 0.7;              % energy conversion coefficient
pN      = 0.1;              % power allocation coefficient
pF      = 1 - pN;
RthN    = 1;                % target data rate of User N bits/s/Hz
RthF    = .1;               % target data rate of User N bits/s/Hz
g2_non  = 2^RthF - 1;
%
SimTimes = 10^4;            % Monte-Carlo repetitions
%
%% Simulation
for ss = 1:length(PS_dB)
    disp(strcat('SNR=',num2str(PS_dB(ss)),'dB'));
    for aa = 1:length(alpha)
        disp(strcat('alpha=',num2str(alpha(aa))));
        for rr = 1:length(rho)
            disp(strcat('rho=',num2str(rho(rr))));
            %
            g2 = 2^(RthF*2/(1-alpha(aa))) - 1; % gamma_2
            
            
            
            
            %alph1 = [0.28195709545229+0.44088121720433i         0.0932454474244731+0.308219428653899i];
             
            alph1 = (0.5)* (rand(1,K)+ 1i*rand(1,K));   % Known transmit correlation coefficent   R_k     
            %Genertion of Co-relation matrix R_k
        RT = [];
         for r1 =1:K
        
                    for ind1 = 1:1:N
                           r(ind1) = alph1(r1)^(ind1-1);
                    end
                    
         R2=toeplitz(r);
         RT = [RT R2];
         clear r r1 ind1 R2
         end         
         
         
         %Transmit beamvectors
    count2 = 1;
    count3 = 1;
    wk=[];
    vk=[];
    for k = 1:K
        RTk = RT(:,count2:count2+(N)-1);
        count2 = count2 + N;
        count3 = count3 + M;
        
        [U1 V1] = eig(RTk);
        wk = [wk U1(:,end)];    %Transmit
        
    end
WK = wk;



count2 = 1;
count3 = 1;
NewTerm=0;
count4=1;
    

    
RTkN = sqrtm(RT(:,count4:count4+(N)-1)); %Picking R_k for a given user 
count4 = count4 + N;
RTkF = sqrtm(RT(:,count4:count4+(N)-1)); %Picking R_k for a given user

         

 for k = 1



        count2 = count2 + N;
        count3 = count3 + M;
        

                
        X1 = WK;
        wk = X1(:,k); %for user desired
        
        wi = X1;   % interference 
        wi(:, k) = [];

        
            A1 =   kron(eye(M), (RTkN*wk).')'*(kron(eye(M), (RTkN*wk).'));%*(1-rho)*pF*PS(ss);
             B2 =  kron(eye(M), (RTkF*wk).')'*(kron(eye(M), (RTkF*wk).'));%*pN*PS(ss);
%             A3=kron(eye(M), (RTkF*wk).')'*(kron(eye(M), (RTkF*wk).'))*(eta.*PS(ss));%*((2*alpha(1)/(1-alpha(1))));
            
%             A1=(RTkN*wk).*(RTkN*wk).';
           
            B1=zeros(M*N,M*N);
            A2=zeros(M*N,M*N);

    for intf = 1:1:size(wi,2) 
            wi2 = (kron(eye(M), (RTkN*wi(:,intf)).')')*(kron(eye(M), (RTkN*wi(:,intf)).'));%*(1-rho)*pN*PS(ss);
            B1 = B1+  wi2;    
            wi3 = (kron(eye(M), (RTkF*wi(:,intf)).')')*(kron(eye(M), (RTkF*wi(:,intf)).'));%*pF*PS(ss);
             A2 = A2+  wi3;

    end
            
            
            
            
            
            
            % channel modelling
            for ii = 1:N
                hSiF(:,ii) = sqrt(lSF/2)*...
                    (randn(SimTimes,1) + 1i*randn(SimTimes,1));
                hSiN(:,ii) = sqrt(lSN/2)*...
                    (randn(SimTimes,1) + 1i*randn(SimTimes,1));
            end
            hNF = sqrt(lNF/2)*...
                (randn(SimTimes,1) + 1i*randn(SimTimes,1));
            
            
              gSiN = hSiN;
            gSiF = hSiF;
            gNF  = abs(hNF.^2);
            
            
            % gsN1=[];
            for i=1:length(gSiN)
               gsN1(i) = (gSiN(i,:)*A1*gSiN(i,:)');
            end
            gSsN=transpose(gsN1);
            
             for i=size(gSiF)
               gSsF1(i) = gSiF(i,:)*B2*gSiF(i,:)'; 
            end
           gSsF = transpose(gSsF1');
            
            % channel gains
            gSiN_ = abs(hSiN.^2);
            gSiF_= abs(hSiF.^2);
            gNF_ = abs(hNF.^2);
            % random selection
            for yy = 1:SimTimes
                i_rand = randperm(K,1);
                gSsN_(yy,1) = gSiN_(yy,i_rand);
                gSsF_(yy,1) = gSiF_(yy,i_rand);
            end
            % Non-cooperative
            snrSsF_non = pF.*PS(ss).*gSsF./(pN.*PS(ss).*gSsF + naF + ncF);   % SINR at F
            % SNR modelling
            snrSsN_xF = (1-rho(rr)).*pF.*PS(ss).*gSsN./...
                ((1-rho(rr)).*pN.*PS(ss).*gSsN ...
                + (1-rho(rr))*naN + ncN);                                           % SINR at N, When msg of user F is detected
            %
            snrSsN_xN = (1-rho(rr)).*pN.*PS(ss).*gSsN/...
                (1-rho(rr))*naN + ncN;                                                  % SINR at N, When msg of User N is detected
            %
            snrSsF = pF.*PS(ss).*gSsF./(pN.*PS(ss).*gSsF + naF + ncF);   % SINR at F for cooperative relaying
            %
            snrNF = eta.*PS(ss).*gSsN.*gNF.*...
                (2*alpha(aa)/(1-alpha(aa))+rho(rr))/(naF + ncF);              % SINR from User N to User F
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%                       Cooperative NOMA
            
            
              snrSsN_xF_ = (1-rho(rr)).*pF.*PS(ss).*gSsN_./...
                ((1-rho(rr)).*pN.*PS(ss).*gSsN_ ...
                + (1-rho(rr))*naN + ncN);                                           % SINR at N, When msg of user F is detected
            %
            snrSsN_xN_ = (1-rho(rr)).*pN.*PS(ss).*gSsN_/...
                (1-rho(rr))*naN + ncN;                                                  % SINR at N, When msg of User N is detected
            %
            snrSsF_ = pF.*PS(ss).*gSsF_./(pN.*PS(ss).*gSsF_ + naF + ncF);   % SINR at F for cooperative relaying
            %
            snrNF_ = eta.*PS(ss).*gSsN_.*gNF_.*...
                (2*alpha(aa)/(1-alpha(aa))+rho(rr))/(naF + ncF);              % SINR from User N to User F
            % count outage events
            count_coop_ = 0;
            
            
            
            
            
            % count outage events
            count_coop = 0;
            count_non = 0;
            %
            for zz = 1:SimTimes
                % for non-cooperative
                if (snrSsF_non(zz) < g2_non)
                    count_non = count_non + 1;
                end
                % for Cooperative 
                if (snrSsN_xF(zz) >= g2) && ...
                        (max(snrSsF(zz),snrNF(zz)) < g2)
                    count_coop = count_coop + 1;
                elseif (snrSsN_xF(zz) < g2) && (snrSsF(zz) < g2)
                    count_coop = count_coop + 1;
                end
            end
            
             for zz = 1:SimTimes
               
              
                % for Cooperative 
                if (snrSsN_xF_(zz) >= g2) && ...
                        (max(snrSsF_(zz),snrNF_(zz)) < g2)
                    count_coop_ = count_coop_ + 1;
                elseif (snrSsN_xF_(zz) < g2) && (snrSsF_(zz) < g2)
                    count_coop_ = count_coop_ + 1;
                end
            end
            
            
            OP_coop_random(aa,rr) = count_coop/SimTimes;
            OP_non_random(aa,rr) = count_non/SimTimes;
            OP_coop_random_(aa,rr) = count_coop_/SimTimes;
        end
    end
    end
end

%% plot

surf(alpha,rho,OP_coop_random,'linestyle',':')
hold on

surf(alpha,rho,OP_coop_random_,'linestyle','-.')

%hold on
%surf(alpha,rho,OP_non_random,'linestyle','-')
%hold on

%
% zlim([10^-3 10^0])
set(gca, 'ZScale', 'log')
%
xlabel('\rho')
ylabel('\alpha')
zlabel('Outage Probability')
legend('Proposed Method','TAS method')%,'Without Cooperatine NOMA.')
%
set(gca,'XTick',0:.2:1)
set(gca,'YTick',0:.2:1)
