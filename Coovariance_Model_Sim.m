%Simulation Setup
function[A, B, Y, CDF] = Coovariance_Model_Sim(wk, M, N, K, var_v, Phi_k,Phi_F,Phi_NF,PS,PS_dB,alpha,rho,RthF,naN ,ncN,pN,pF,ncF,naF,eta, NoumberOfTrial,xbins)
Y = [];
A = [];
B = [];
CDF = [];
count2 = 1;
count3 = 1;
    
    for k = 1:K
        

%                 
        X1 = wk;
        wk_k = X1(:,k);
        
        wi = X1;
        wi(:, k) = [];
        
           vec1 = sqrtm(Phi_k{k})*wk_k;
           vec1_1 = sqrtm(Phi_F{k})*wk_k;
           vec1_2 = sqrtm(Phi_NF{k})*wk_k;
           
           A1 =   vec1*vec1';
           A1_1 =   vec1_1*vec1_1';
           A1_2 =   vec1_2*vec1_2';
%       
%             keyboard
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
            Bk=zeros(M,M);
             Bk_1=zeros(M,M);
              Bk_2=zeros(M,M);

            for intf = 1:size(wi,2) 
            
                vec2 = sqrtm(Phi_k{k})*wi(:,intf);                             %sqrt(1/(M*var_v))*(kron(eye(M), (RTk*wi(:,intf)).')'*vk_B')' ;
                vec2_1 = sqrtm(Phi_F{k})*wi(:,intf); 
                vec2_2 = sqrtm(Phi_NF{k})*wi(:,intf); 
                
                wi2  = vec2*vec2';
                wi2_1  = vec2_1*vec2_1';
                wi2_2  = vec2_2*vec2_2';
                
       
            Bk = Bk+  wi2;
            Bk_1 = Bk_1+  wi2_1;
            Bk_2 = Bk_2+  wi2_2;
            end
            B1 = Bk;
            B1_1 = Bk_1;
            B1_2 = Bk_2;
%        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                          





for ss = 1:length(PS_dB)
    count_2 = 0;
    disp(strcat('SNR=',num2str(PS_dB(ss)),'dB'));
    for aa = 1:length(alpha)
        disp(strcat('alpha=',num2str(alpha(aa))));
        for rr = 1:length(rho)
            disp(strcat('rho=',num2str(rho(rr))));
            %
            g2 = 2^(RthF*2/(1-alpha(aa))) - 1; % gamma_2
            
%%%%%%%%%%%%%%%%%%%%
%%%   Channels   %%%

h_N1=sqrt(1/2).*(randn(NoumberOfTrial,M)+ 1i*randn(NoumberOfTrial,M));     %Row vector of length N
h_F1=sqrt(1/2)'.*(randn(NoumberOfTrial,M)+ 1i*randn(NoumberOfTrial,M));     %Row vector of length N
h_NF1=sqrt(1/2)'.*(randn(NoumberOfTrial,M)+ 1i*randn(NoumberOfTrial,M));     %Row vector of length N
hhN=sqrt(1/2)'.*(randn(NoumberOfTrial,M)+ 1i*randn(NoumberOfTrial,M));
hhF=sqrt(1/2)'.*(randn(NoumberOfTrial,M)+ 1i*randn(NoumberOfTrial,M));
hhNF=sqrt(1/2)'.*(randn(NoumberOfTrial,M)+ 1i*randn(NoumberOfTrial,M));
   
 
      for runs=1:NoumberOfTrial
          h_N=h_N1(runs,:);
                    %h_N=transpose(h_N_');

          h_F=h_F1(runs,:);
                    %h_F=transpose(h_F_');

         h_NF=h_NF1(runs,:);
         %h_NF_= (h_NF)
          gNF= abs(h_NF.^2);
        
        %Z_N=(rBF*hhN(runs,:)*diag(Lamc)*hhN(runs,:)');
        
%        Z_N=rBF*hhN(runs,:)*diag(Lamc)*(rBN*hhN(runs,:))'; 
%        Z_F=rBF*hhF(runs,:)*diag(Lamd)*(rBF*hhN(runs,:))';
%        Z_NF=rBF*hhNF(runs,:)*diag(Lamc)*(rNF*hhN(runs,:))'; 


       
       SINR_N_sF = (1-rho(rr)).*pF.*PS(ss).*(h_N*A1*h_N')./ ((1-rho(rr)).*pN.*PS(ss).*real(h_N*B1*h_N') + (1-rho(rr))*naN + ncN);
        SINR_F_sF =pF.*PS(ss).*(h_F*A1*h_F')./(pN.*PS(ss).*(h_F*B1*h_F') + naF + ncF );
        SINR_NF = eta.*PS(ss).*(h_N*A1*h_N').*h_NF.* (2*alpha(aa)/(1-alpha(aa))+rho(rr))/(naF + ncF);
   
 
 if SINR_N_sF >= g2 && SINR_F_sF< g2
     count_2 = count_2 + 1;
     
 elseif (SINR_N_sF < g2) && (SINR_F_sF < g2)
     count_2 = count_2 + 1;
 end
      end
      CDF(ss) = count_2./NoumberOfTrial;
    %  OP_S2_F_sim1(ss) = count_2./NoumberOfTrial;
        
 %   [envlp,CDF]=hist(real(OP_S2_F_sim1),100);
            
 
            
            
            
            
      end
        end
end





























% for runs=1:NoumberOfTrial
% 
%                                h=sqrt(1/2)'.*(randn(M,1)+ 1i*randn(M,1));     %Row vector of length N
%                                Y1(runs)=real((((h'*A1*h)/((var_v)+h'*B1*h))));
% 
%                             end
%                    
%                     [envlp,y]=hist(real(Y1),xbins);
%                     W=y(2)-y(1);
%                     pdf_Y=envlp/(length(Y1)*W);
%                     CDF1 =(cumsum(W*pdf_Y)); 
%                      
%                     
%                             A = [A A1];
%                             B = [B B1];
%                             Y = [Y Y1'];
%                             CDF = [CDF CDF1'];
%                     clear W   
    end
end