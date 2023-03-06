function[CDFA] = Coovariance_Model_Analy(wk, M, N, K, var_v, Phi_k,PS,PS_dB,alpha,rho,RthF,pF,pN, tol)

A = [];
B = [];
 
    for k = 1:K
        

%                 
        X1 = wk;
        wk_k = X1(:,k);
        
        wi = X1;
        wi(:, k) = [];
        
           vec1 = sqrtm(Phi_k{k})*wk_k;
           A1 =   vec1*vec1';

            Bk=zeros(M,M);

            for intf = 1:size(wi,2) 
            
                vec2 = sqrtm(Phi_k{k})*wi(:,intf);                             %sqrt(1/(M*var_v))*(kron(eye(M), (RTk*wi(:,intf)).')'*vk_B')' ;
                wi2  = vec2*vec2';
       
            Bk = Bk+  wi2;        
            end
            B1 = Bk;
                            A = [A A1];
                            B = [B B1];
    end


T = M;
CDFA = [];

    count2 = 1;
    for k = 1:1
        
        A1A = A(:,count2:count2+(M)-1);
        B1A = B(:,count2:count2+(M)-1);
        count2 = count2 + (M);
%         keyboard
        
        
           %     for kk=1:length(y)
 
           
           for ss = 1:length(PS_dB)
    count_2 = 0;
    disp(strcat('SNR=',num2str(PS_dB(ss)),'dB'));
    for aa = 1:length(alpha)
        disp(strcat('alpha=',num2str(alpha(aa))));
        for rr = 1:length(rho)
            disp(strcat('rho=',num2str(rho(rr))));
            %
            g2 = 2^(RthF*2/(1-alpha(aa))) - 1; % gamma_2
           
           
            dSN     = 3;
            dSF     = 10; 
            dNF     = dSF - dSN;
            L       = 1e3;  
            lSN     = L*dSN^-3;  
            lNF     = L*dNF^-3;
            mu_a = g2/(1-rho(rr))/(pF-pN*g2);
            eta     = 0.8; 
            naN     = (10^(-7))*1e6;    % naN = -100 dBm, BW = 1 MHZ_N
ncN     = (10^(-6))*1e6;    % naN = -90 dBm,  BW = 1 MHZ_N
naF     = (10^(-7))*1e6;
ncF     = (10^(-6))*1e6;
             c  = eta*PS(ss)*(2*alpha(aa)/(1-alpha(aa))+rho(rr))/(naF + ncF);
            
                        P= A1.*pF.*PS(ss)  - g2*B1.*pN.*PS(ss);      %(2^(y(kk))-1)
                         P1= B1.*pF.*PS(ss)  - g2*A1.*pN.*PS(ss);      %(2^(y(kk))-1)
                        Pmat{k} =P;
%                         [Up,Lamp]=eig(P);
% 
%                         lamp=(diag(Lamp));
%                         lamp = sort(lamp,'descend')';
% 
%                         lamp2 = lamp;
%                         L1 = size(lamp(lamp >tol));
%                         T1 = size(lamp,2);
%                         for n = 1:L1(2)
%                             lamp1 = lamp2;
%                             lamp1(n)= [];                             
%                             for j = 1: T1-1                        
%                                 Prod1(j) = lamp(n)- lamp1(j);                            
%                             end                            
%                             Prod    = prod(Prod1);                            
%                             Sum1(n) = ((lamp(n))^(T) /abs(lamp(n))*(1/Prod(n))*exp(-g2*(var_v)/lamp(n))*heaviside(g2*(var_v)/lamp(n)));
%                         end
%                         
%                       Sum = sum(Sum1);
%                       
%                       FA(ss)= heaviside(g2*(var_v))-Sum;     
%                       
%                     extpara = lSN.*trace(A1A);
%                     expterm(ss) =exp(-mu_a/extpara); 
%                      int(ss) = integral(@(t) exp(-t-g2/(extpara*lNF*c*t)),mu_a/extpara,Inf,'ArrayValued',true);
%                      newterm(ss)= expterm(ss) - int(ss) ;
%                       FA3(ss) = FA(ss)*newterm(ss);
%                      output(ss) = FA(ss) .*(FA_1(ss)+ newterm(ss));  
%                       
%                       
%                       
%                       
%                       CDFA_1 = [CDFA FA3(ss)];
%                       CDFA = [CDFA FA(ss)];
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
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
                        
                        T1 = size(lamp,2);
                        for n = 1:L1(1)
                        
                            lamp1 = lamp2;
                             lamp11 = lamp_2;
                            
                            lamp1(n)= [];
                            lamp11(n)= [];
                            
                            for j = 1: T1-1;                        
                                Prod1(j) = lamp(n)- lamp1(j);    
                                Prod11(j) = lamp_1(n)- lamp11(j);  

                                
                            end                            
                            Prod = prod(Prod1);     
                            Prod_1 = prod(Prod11); 
                            
                      Sum1(n) = ((lamp(n))^T /abs(lamp(n))*(1/Prod(n))*exp(-g2*(naF+ncF)/lamp(n))*heaviside(g2*(naF+ncF)/lamp(n)));
                      Sum2(n) = ((lamp_1(n))^T /abs(lamp_1(n))*(1/Prod_1(n))*exp(-g2*((1-rho).*naN+ncN)/lamp_1(n))*heaviside(g2*((1-rho).*naN+ncN)/lamp_1(n)));

                                
                                
                        end
                        
                      Sum = sum(Sum1);
                      Sum_1 = sum(Sum2);


                     extpara = lSN.*trace(A1);


                      FA(ss)= heaviside(g2*(naF+ncF))-Sum;  

                      FA_1(ss)= heaviside(g2*((1-rho).*naN+ncN))-Sum_1;  

                      
                      
                        expterm(ss) =exp(-mu_a/extpara); 
                        
                       int(ss) = integral(@(t) exp(-t-g2/(extpara*lNF*c*t)),mu_a/extpara,Inf,'ArrayValued',true)

                       newterm(ss)= expterm(ss) - int(ss) ;

                       FA3(ss) = FA(ss)*newterm(ss);

                       output(ss) = FA(ss) .*(FA_1(ss)+ newterm(ss));               %Outage Probability of user-B
                       FA_(ss) = output(ss);                                     %Outage Probability of user-A
                       
                      
                       CDFA = [CDFA FA(ss)];
                      
                      
                     
                end
                
                
       
    end
    
           end
           
    end
    
end
