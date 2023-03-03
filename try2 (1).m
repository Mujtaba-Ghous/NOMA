%close all
%clear all

hold on
%
%% Simulation parameters
%
N       = 7;                      % # of antenna
M = 1; % Number of recieve antennas
K = 2; % Numer of users at transmit side
rho     = .3;         % power splitting ratio
alpha   = 0.5;         % time fraction for EH
PS_dB   = -30:2:20;             % transmit SNR = Ps/N0 in dB
PS      = 10.^(PS_dB./10);
naN     = (10^(-7))*1e6;    % naN = -100 dBm, BW = 1 MHz
ncN     = (10^(-6))*1e6;    % naN = -90 dBm,  BW = 1 MHz
naF     = (10^(-7))*1e6;
ncF     = (10^(-6))*1e6;
epsilon = 3;                % pathloss exponent
dSF     = 10;                   % S-F distance
dSN     = 3;
dNF     = dSF - dSN;
L       = 1e3;                  % path-loss at reference distance
%


lSN     = L*dSN^-3;         %1.5%    % lambda
lSF     = L*dSF^-3;
lNF     = L*dNF^-3;
%
eta     = 0;%0.7;              % energy conversion coefficient
RthN    = .17;                % target data rate of User N bits/s/Hz
RthF    = .17;               % target data rate of User N bits/s/Hz
%[pN,pF] = PowerAllocation(RthN,RthF);
%
pN = (2^(2*RthN) - 1)/(2^(2*RthN+2*RthF)-1);
pF = 1-pN;
SimTimes = 1000;           % Monte-Carlo repetitions
%
 rng(1,'twister') 
 tol  = 1e-12
%% Simulation
%
for ss = 1:length(PS_dB)
    disp(strcat('SNR=',num2str(PS_dB(ss)),'dB'));
    for aa = 1:length(alpha)
        disp(strcat('alpha=',num2str(alpha(aa))));
        for rr = 1:length(rho)
            disp(strcat('rho=',num2str(rho(rr))));
            %
            g2 = 2^(RthF*2/(1-alpha(aa))) - 1; % gamma_2
        
          
            %alph1 = [7.78195709545229+0.44088121720433i         5.1932454474244731+0.308219428653899i];
            
        % alph1 = [0.25195709545229+0.44088121720433i         0.0932454474244731+0.308219428653899i]; % for k=7
        alph1 = [0.28195709545229+0.44088121720433i  0.0932454474244731+0.308219428653899i]; % for k=5
         
         %   alph1 = [0.20195709545229+0.44088121720433i         0.0932454474244731+0.308219428653899i]; % for k=11
           
           %alph1 = [0.37195709545229+0.44088121720433i         0.0332454474244731+0.308219428653899i];  %for K=3 
              
              
           
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
                hiF(:,ii) = sqrt(lSF)*...
                    (randn(SimTimes,1) + 1i*randn(SimTimes,1));
                hiN(:,ii) = sqrt(lSN)*...
                    (randn(SimTimes,1) + 1i*randn(SimTimes,1));
            end
            H1=hiF.';
            H1_=hiF;
            H2=hiN.';
            H2_=hiN;
           
           
            hNF = sqrt(lNF/2)*...
                (randn(SimTimes,1) + 1i*randn(SimTimes,1));
             H3=hNF.';
            H3_=hNF;
%             H11=H1_*H1;
%             H12=inv(H11);
%             H13=H1*H12;
            H14=H1*inv(H1_*H1);
            H15=transpose(H14);
           %hiF=H15;
            
            H16=H2*inv(H2_*H2);
            H17=transpose(H16);
            %hiN=H17;
            
            H18=H3*inv(H3_*H3);
            H19=transpose(H18);
          %hNF= H19;
            
            %H11=H1.*inv(H1_.*H1);
            %H22=
            %H33=
            
            % channel gains
            giN = (hiN);
            giF = (hiF);
            gNF = abs(hNF.^2);
            A1  = (A1);
            
           % gsN1=[];
            for i=1:length(giN)
               gsN1(i) = (giN(i,:)*A1*giN(i,:)');
            end
            gsN=transpose(gsN1);
            
             for i=size(giF)
               gSsF1(i) = giF(i,:)*B2*giF(i,:)'; 
            end
           gSsF = transpose(gSsF1');
            
%             find the best far
%          [gSsF,I] = max(giF,[],2);
%             for yy = 1:SimTimes
%                 gsN(yy,1) = giN(yy,I(yy));
%             end
%            SNR modelling
             snrsNxF = (1-rho(rr)).*pF.*PS(ss).*gsN./ ((1-rho(rr)).*pN.*PS(ss).*gsN + (1-rho(rr))*naN + ncN);
%             %
             snrsNxN = (1-rho(rr)).*pN.*PS(ss).*gsN/(1-rho(rr))*naN + ncN;
%             %
         snrsF = pF.*PS(ss).*gSsF./(pN.*PS(ss).*gSsF + naF + ncF);
%             %
            snrNF = eta.*PS(ss).*gsN.*gNF.* (2*alpha(aa)/(1-alpha(aa))+rho(rr))/(naF + ncF);
%             % count outage events
             count_2 = 0;

            for zz = 1:SimTimes
                if (snrsNxF(zz) >= g2) && ...
                (snrsF(zz)> g2)
                %   (max(snrsF(zz),snrNF(zz)) < g2)   
              
                    count_2 = count_2 + 1;
                elseif (snrsNxF(zz) < g2) && (snrsF(zz) < g2)
                    count_2 = count_2 + 1;
                end
            end
            OP_S2_F_sim1(ss) = count_2./SimTimes;
        
            
            %% Analytical Results
            a1 = (1-rho(rr))*pF*PS(ss)/((1-rho(rr))*naN + ncN);
            a2 = (1-rho(rr))*pN*PS(ss)/((1-rho(rr))*naN + ncN);
            b1 = pF * PS(ss) / (naF + ncF);
            b2 = pN * PS(ss) / (naF + ncF);
            c  = eta*PS(ss)*(2*alpha(aa)/(1-alpha(aa))+rho(rr))/(naF + ncF);
            %
            mu_a = g2/(a1-a2*g2);
            mu_b = g2/(b1-b2*g2);
            %
            c_ = 0*(2*alpha(aa)/(1-alpha(aa))+rho(rr));
            mu_a_ = g2/(1-rho(rr))/(pF-pN*g2);
            SNR = PS(ss)/(naN+ncN);
            
            
                          P= A1.*(1-rho(rr)).*pF.*PS(ss) - g2.*A1*(1-rho(rr)).*pN.*PS(ss);      
                          P1= B2.*pF.*PS(ss) - g2.*B2*pN.*PS(ss);      

                          T = M*N;
                          Prod1=[];
                          Prod11=[];
                          

                        [Up,Lamp]=eig(P);
                        [Up,Lamp_1]=eig(P1*lSN);
                        
                        lamp=(diag(Lamp));
                        lamp_1=(diag(Lamp_1));
                        
                        lamp = sort(lamp,'descend')';
                        lamp_1 = sort(lamp_1,'descend')';
                        
%                         keyboard
                        lamp2 = lamp;
                          lamp_2 = lamp_1;
                        
                        L1 = size(lamp(lamp >tol));
                        L_1 = size(lamp(lamp_1 >tol));
                        
                        T1 = size(lamp,2);
                        T_1 = size(lamp_1,2);
                        for n = 1:L1(1)
                        
                            lamp1 = lamp2;
                             lamp11 = lamp_2;
                            
                            lamp1(n)= [];
                            lamp11(n)= [];
                            
                            for j = 1: T-1;                        
                                Prod1(j) = lamp(n)- lamp1(j);    
                                Prod11(j) = lamp_1(n)- lamp11(j);                            

                            end                            
                            Prod = prod(Prod1);     
                            Prod_1 = prod(Prod11);     
                          %  Sum1(n) = ((lamp(n))^T /abs(lamp(n))*(1/Prod(n))*exp(-y(kk)*(M*var_v)/lamp(n))*heaviside(y(kk)*(M*var_v)/lamp(n)));
                               Sum1(n) = ((lamp(n))^T /abs(lamp(n))*(1/Prod(n))*exp(-g2*(naF+ncF)/lamp(n))*heaviside(g2*(naF+ncF)/lamp(n)));
                                Sum2(n) = ((lamp_1(n))^T /abs(lamp_1(n))*(1/Prod_1(n))*exp(-g2*((1-rho).*naN+ncN)/lamp_1(n))*heaviside(g2*((1-rho).*naN+ncN)/lamp_1(n)));
%                             clear Prod Prod1
                        end
                        
                      Sum = sum(Sum1);
                      Sum_1 = sum(Sum2);

                    extpara = lSN.*trace(A1);
                    % extpara = trace(A1).*lSN;% lSN ;%/ trace(B2); %;.* trace(A1);

                     
                      FA(ss)= heaviside(g2*(naF+ncF))-Sum;    
                        

                     FA_1(ss)= heaviside(g2*((1-rho).*naN+ncN))-Sum_1;  
                      
                      
                        expterm(ss) =exp(-mu_a/extpara); 
                       int(ss) = integral(@(t) exp(-t-g2/(extpara*lNF*c*t)),mu_a/extpara,Inf,'ArrayValued',true)
                       newterm(ss)= expterm(ss) - int(ss) ;
                      % output(ss)=  FA(ss).*(FA_1(ss)+newterm(ss));
                       FA3(ss) = FA(ss)*newterm(ss);
                       output(ss) = FA(ss) .*(FA_1(ss)+ newterm(ss));
                       

                end
      
        end
    end
    end


%% plot
%plot(PS_dB,OP_S2_F_sim1,'o', PS_dB,output,'-');
plot(PS_dB,output,'-');
%plot(PS_dB,OP_S2_F_sim1,'-');


%xlabel('SNR (dBm)')
%ylabel('Outage Probability')
%plot(y,CDFA(:,K),'r');
%plot(PS_dB,output,'o')
%, ...
   % PS_dB,OP_S2_F_ana,'-')
% axis([-20 40 0 1])
%legend('Simulation','Analytical')
%legend('  N=7 (sim.)  ',' N=7 (ana.)  ')



legend('  N=5 Proposed Method (sim.)  ',' N=5 Proposed Method (ana.)','  N=5 TAS Technique [30](sim.)  ',' N=5 TAS Technique [30](ana.)','  N=5 Proposed Method (sim.)  ',' N=5 Proposed Method (ana.)','  N=5 TAS Technique [30](sim.)  ',' N=5 TAS Technique [30](ana.)' )



