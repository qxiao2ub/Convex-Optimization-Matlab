% 0515-Make all var LN and continuous
% 0508-denotes are from E:\OneDrive - University at Buffalo\2-Buffalo\1-Course\1-Kang Research\5-Send to Kang\wk13-0421-meeting\Cont_Extended_Abstract_0424
% 0424-operational formulation explain look at E:\OneDrive - University at Buffalo\2-Buffalo\1-Course\1-Kang Research\3.Zhiheng Data\2020.04.05 Formulation and Code "Formulation.ppt"
% 0508-hide c17 that limit Uijt in operational, otherwise no solution.[this action is correct, this c should put in upper for next stage]

%0-read all test case 12 [done by wk2-02/07/20fri]

%0-start recording running time
tic

%0.1-read and store test case 12
clear;clc;close all;
d= zeros(20,20,40);%create 20*20*40 matrix
% 101520-read different cijt
c= zeros(20,20,40);
% 101520-pricing range for every i,j,t
P= zeros(20,20,40);


% read dijt
for t=1:40 %40 time periods
    fileName=['Period',num2str(t-1),'.txt']; %txt: Period0 to 39
    read=load(fileName); %load that txt
    d(:,:,t)=read; %assign txt to dijt
end

% 101520-read cijt and keep standard $26/hrs
for t1=1:40 %40 time periods
    fileName1=['cijtPeriod',num2str(t1-1),'.txt']; 
    read1=load(fileName1); %load that txt
    c(:,:,t1)=read1./6; %assign txt to dijt
end


% %{
% 0520-test whether z_op decreases as Pijt increases
% for Pijt=1:1:1
%   z_op(Pijt) = zeros(Pijt);
% end

%1-decision variable U111 to U(20,20,40),X111 to X(20,20,40),V1,1 to V(20,40)

%-----for small range before whole dataset, set i,j,t upper bound------
i_up=2; %should be 20
j_up=2; %should be 20
t_up=4; %should be 40
%--------
% 0517-M is c16 time window, set at 1 gives it relaxation
% 0424-M is c16 time window, set at 3 gives it relaxation
M=1; 

%0515-change all var to continours instead of MIP
% X = sdpvar(i_up,j_up,t_up+M); %0424-Xi,j,t+M, cuz time window constraint c16 Xijt is beyond t
% %X = intvar(i_up,j_up,t_up+M); %0424-Xi,j,t+M, cuz time window constraint c16 Xijt is beyond t
% X_t0 = sdpvar(i_up,j_up,1); %0424-Xi,j,0 means moves initially for c18: inventory balance eq
% %X_t0 = intvar(i_up,j_up,1); %0424-Xi,j,0 means moves initially for c18: inventory balance eq
% U = sdpvar(i_up,j_up,t_up);
% %U = intvar(i_up,j_up,t_up); 
% U_t0 = sdpvar(i_up,j_up,1); %0424-this is for Ui,j,0 
% %U_t0 = intvar(i_up,j_up,1); %0424-this is for Ui,j,0 
% V = sdpvar(i_up,t_up); 
% %V = intvar(i_up,t_up); 
% V_t0 = sdpvar(i_up,1); %0424-this is especially for Vi,0
% %V_t0 = intvar(i_up,1); %0424-this is especially for Vi,0
% D_ij_SAV = sdpvar(i_up,j_up,t_up); %0425-this is to show D_ij_SAV, for futher plot, then find upper bound of Pijt
% %D_ij_SAV = intvar(i_up,j_up,t_up); %0425-this is to show D_ij_SAV, for futher plot, then find upper bound of Pijt

%1.1-parameters
% a=b=m=n=1,T=40
% a=1;
% b=1;
% m=1;
% n=1;
w_c=1.3;
w_w=39;
w_p=1.3;
q=2.62;
k=0.053;
T=1; %total T is 40(look at obj func), this relates to t_up
T_abs=1; %temporarily set |T|=2 or 3 is related to game theory that |T|
% cijt=40, pijt=10, hit=10, alphaijt=0.1 (all c,p,h,alpha are subjectively to be same)
% cijt=10;
pijt=1.875; % 092920-this is $1.875/time period = $5/hrs
hit=1;
alpha=1;

%initially price Pijt is very low, set to 1
%Pijt=1; %0520 - temporarily hide to see what happens as Pijt increases for oper_obj

% fid = fopen(['Pijt1to10001by1000PerStep_OperCostSplit.txt'],'wt'); %creat & write z_op to txt

%1.2-add operational level constraints

% Pijt=13;

% c1=0;

% CompU=zeros(16000,2);

% P_ij_SAV= zeros(i_up,j_up,t_up);

% for P_pa=1.3:0.1:1.4
%   P_pa 
% end


% %{

for P_pa=1:1 %100820-$6 to $24/TP = $17 to $68/hrs

% Pijt=13;

% c1=c1+1;

X = sdpvar(i_up,j_up,t_up+M); %0424-Xi,j,t+M, cuz time window constraint c16 Xijt is beyond t

X_t0 = sdpvar(i_up,j_up,1); %0424-Xi,j,0 means moves initially for c18: inventory balance eq

U = sdpvar(i_up,j_up,t_up);

U_t0 = sdpvar(i_up,j_up,1); %0424-this is for Ui,j,0 

V = sdpvar(i_up,t_up); 

V_t0 = sdpvar(i_up,1); %0424-this is especially for Vi,0

D_ij_SAV = sdpvar(i_up,j_up,t_up); %0425-this is to show D_ij_SAV, for futher plot, then find upper bound of Pijt

% 101520-adjust Pijt based on ratio
for i=1:i_up
  for j=1:j_up
    for t=1:t_up
      if c(i,j,t)<=5
        P(i,j,t)=7;
      elseif c(i,j,t)>5 && c(i,j,t)<=10
        P(i,j,t)=14;
      elseif c(i,j,t)>10 && c(i,j,t)<=15
        P(i,j,t)=21;
      else 
        P(i,j,t)=28;
      end
    end
  end
end

% 100820-express D_ij_SAV
for i=1:i_up
  for j=1:j_up
    for t=1:t_up
      D_ij_SAV(i,j,t)=((k.*w_p.*P(i,j,t)+c(i,j,t).*w_c.*q).*T_abs.*U(i,j,t))./(c(i,j,t).*d(i,j,t).*w_c.*q-k.*w_w.*T_abs.*U(i,j,t));
%       D_ij_SAV(i,j,t)=((k.*w_p.*Pijt+cijt.*w_c.*q).*T_abs.*U(i,j,t))./(cijt.*d(i,j,t).*w_c.*q-k.*w_w.*T_abs.*U(i,j,t));
    end
  end
end

% end



% %{

C = [];

% c26-092120overleaf
for i=1:i_up %c26
   for j=1:j_up
     for t=1:t_up
       if (t==1)
%          C = [C, U(i,j,1) >= U_t0(i,j,1)+(b*d(i,j,1)+cijt*d(i,j,1)*m-n*T_abs*U(i,j,1))./(a+b+m*(cijt+Pijt))-X(i,j,1)];
         C = [C, U(i,j,t) >= U_t0(i,j,t)+(c(i,j,t)*d(i,j,t)*w_c*q-k*w_w*T_abs*U(i,j,t))./(k*w_p*P(i,j,t)+c(i,j,t)*w_c*q)-X(i,j,t)]; %c26: t=1
       else       
         C = [C, U(i,j,t) >= U(i,j,t-1)+(c(i,j,t)*d(i,j,t)*w_c*q-k*w_w*T_abs*U(i,j,t))./(k*w_p*P(i,j,t)+c(i,j,t)*w_c*q)-X(i,j,t)]; %c26: t from 2
       end
     end
   end    
end

% 0424-c16-I suppose time window is M period for c16 - relax this one, hope get solution
% this also considers Xijt beyond T
for i=1:i_up %c27
    for j=1:j_up
      for t=1:t_up
        if (t==1)
%           C = [C, sum(X(i,j,t:(t+M))) >= U_t0(i,j,1)+(b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t))./(a+b+m*(cijt+Pijt))]; %c16: t=1
          C = [C, sum(X(i,j,t:(t+M))) >= U_t0(i,j,t)+(c(i,j,t)*d(i,j,t)*w_c*q-k*w_w*T_abs*U(i,j,t))./(k*w_p*P(i,j,t)+c(i,j,t)*w_c*q)]; 
        else
          C = [C, sum(X(i,j,t:(t+M))) >= U(i,j,t-1)+(c(i,j,t)*d(i,j,t)*w_c*q-k*w_w*T_abs*U(i,j,t))./(k*w_p*P(i,j,t)+c(i,j,t)*w_c*q)]; 
        end
      end  
    end     
end

% 0424 - c17 is in upper level, no need in lower level
% for i=1:i_up %c17
%   for j=1:j_up
%     for t=1:t_up
%       C = [C, ((b-1)*d(i,j,t)+cijt*d(i,j,t)*m)./(n*T_abs) <= U(i,j,t)];
%       C = [C, U(i,j,t) <= (b*d(i,j,t)+cijt*d(i,j,t)*m)./(n*T_abs)];
%     end
%   end
% end

% 0424-c18-here makes <= for Vi,t cuz alpha all set to 1, too relaxed, should be variaty
for i=1:i_up %c28
  for t=1:t_up
    if (t==1)
      C = [C, V(i,t) == V_t0(i,t) + sum(alpha*X_t0(1:j_up,i,t)) - sum(X(i,1:j_up,t)) ]; 
    else
      C = [C, V(i,t) == V(i,t-1) + sum(sum(alpha*X(1:j_up,i,1:(t-1)))) - sum(X(i,1:j_up,t)) ]; %(remember to use 2 sum for 2nd)
    end
  end
end

% 0508-c29-not added, cuz it looks already relaxed in math function, hidden property of our formula

for i=1:i_up %c30: Xijt
    for j=1:j_up
        for t=1:(t_up+M)
          C = [C, X(i,j,t) >= 0];
        end
    end    
end    

for i=1:i_up %c20: X_t0(i,j,1) is Xijt=0
    for j=1:j_up
        C = [C, X_t0(i,j,1) >= 0];
    end    
end    

for i=1:i_up %c20: Uijt
    for j=1:j_up
        for t=1:t_up
            C = [C, U(i,j,t) >= 0];
        end
    end    
end    

for i=1:i_up %c20: U_t0(i,j,1) is Uijt=0
    for j=1:j_up
        C = [C, U_t0(i,j,1) >= 0];
    end    
end 

for i=1:i_up %c20: Vit
    for t=1:t_up
      C = [C, V(i,t) >= 0];
    end    
end    

for i=1:i_up %c20: V_t0(i,1) is V(i,0)
  C = [C, V_t0(i,1) >= 0];
end

% 092220-Uijt<=dijt
for i=1:i_up %c20: Uijt
    for j=1:j_up
        for t=1:t_up
            C = [C, U(i,j,t) <= d(i,j,t)];
        end
    end    
end    

% 100820-Uijt<=Dijt
% for i=1:i_up 
%     for j=1:j_up
%         for t=1:t_up
%             C = [C, U(i,j,t) <= D_ij_SAV(i,j,t)];
%         end
%     end    
% end    

% 092920-Xijt<=dijt
% for i=1:i_up %c20: Uijt
%     for j=1:j_up
%         for t=1:t_up
%             C = [C, X(i,j,t) <= d(i,j,t)];
%         end
%     end    
% end    

%1.3-solve operational level: settings: use cplex as solver

%0515-following is set optimality for MIP, but kang said make it LNContinuous-https://www.ibm.com/support/knowledgecenter/SSSA5P_12.9.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/OptimalityTarget.html
%0515-'cplex.optimalitytarget'=1: Searches for a globally optimal solution to a convex model.
ops = sdpsettings ('solver','cplex','verbose',2,'cplex.optimalitytarget',1);%0508-verbose 2 shows iterations, 1 means shown normal message, 0 shows silent message
%0515-MIP limit tolerance gap to 1-https://www.ibm.com/support/knowledgecenter/SSSA5P_12.9.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/EpGap.html
%ops = sdpsettings ('solver','cplex','verbose',2,'cplex.mip.tolerances.mipgap',1);%0508-verbose 2 shows iterations, 1 means shown normal message, 0 shows silent message

%0306-limit barrier iteration-https://www.ibm.com/support/knowledgecenter/SSSA5P_12.9.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/BarItLim.html
%Cplex.Param.barrier.limits.iteration = 4;[not work]

%0306-limit network simplex iteration-https://www.ibm.com/support/knowledgecenter/SSSA5P_12.9.0/ilog.odms.cplex.help/CPLEX/Parameters/topics/NetItLim.html
%Cplex.Param.network.iterations = 2;[not work]

%ops.mip.limits.nodes = 4;[not work]

%0304-use the online method for out-of-memory: https://www.ilovematlab.cn/thread-266244-1-1.html
%ops.mip.strategy.file = 3; %[not work: should put in sdpsetting] this has same results 1159110 node, but reduces time for 3 min.

%0305-use online method for out-of-memory: https://www.ibm.com/developerworks/community/forums/html/topic?id=e3ea344c-7249-4edd-a47d-29e1f80fa480
%Cplex.Param.mip.limits.auxrootthreads = 2; %[not work: should put in sdpsetting] -1 or 1

%obj-operational obj z_op
z_op = sum(sum(sum(c(i,j,t).*X(1:i_up,1:j_up,1:t_up+M)))) + sum(sum(c(i,j,t).*X_t0(1:i_up,1:j_up,1)))...
  + sum(sum(sum(pijt.*U(1:i_up,1:j_up,1:t_up)))) + sum(sum(pijt.*U_t0(1:i_up,1:j_up,1)))...
  + sum(sum(hit.*V(1:i_up,1:t_up))) + sum(hit.*V_t0(1:i_up,1));


% z_op(Pijt) = sum(sum(sum(c(i,j,t).*X(1:i_up,1:j_up,1:t_up+M)))) + sum(sum(c(i,j,t).*X_t0(1:i_up,1:j_up,1)))...
%   + sum(sum(sum(pijt.*U(1:i_up,1:j_up,1:t_up)))) + sum(sum(pijt.*U_t0(1:i_up,1:j_up,1)))...
%   + sum(sum(hit.*V(1:i_up,1:t_up))) + sum(hit.*V_t0(1:i_up,1));

%1.35-express D_ij_SAV
%   for i=1:i_up
%     for j=1:j_up
%       for t=1:t_up
%         D_ij_SAV(i,j,t) = (b*d(i,j,t)+cijt*d(i,j,t)*m-n*T_abs*U(i,j,t))./(a+b+m*(cijt+Pijt));
%       end
%     end
%   end

%1.4-check if solved and output
result = optimize(C,z_op,ops); % 101520
% result = optimize(C,z_op(Pijt),ops); % 0508-origin
% result = solvesdp(C,z_op,ops); % 0508-this works

%   result = optimize(C,z);
if result.problem == 0 % problem =0 means solve successfully,only print Uijt here

    fprintf('\nz_op(Pijt) value: ')
    value(z_op)
%     value(z_op(Pijt))
    
    fprintf('\nPijt value: ')
    value(P_pa)
    
%     fprintf(fid,'%.4f\t',double(z_op(Pijt))); %0523-save total cost;keep 4 decimals
%     fprintf(fid,'%.4f\t',double(sum(sum(sum(cijt.*X(1:i_up,1:j_up,1:t_up)))) + sum(sum(cijt.*X_t0(1:i_up,1:j_up,1))))); 
%     fprintf(fid,'%.4f\t',double(sum(sum(sum(pijt.*U(1:i_up,1:j_up,1:t_up)))) + sum(sum(pijt.*U_t0(1:i_up,1:j_up,1)))));
%     fprintf(fid,'%.4f\n',double(sum(sum(hit.*V(1:i_up,1:t_up))) + sum(hit.*V_t0(1:i_up,1))));
    
%   %     fprintf('\nD_ij_SAV value: ') %0515-as for now hide D_ij_SAV
%   %     value(D_ij_SAV) 
%   %  
%     fprintf('\nU value: ')
%     value(U)
    
    % 0517-following is save Uijt in txt
%1.6-save Uijt as txt for upper use
% for U_t=1:t_up 
%      fid1 = fopen(['Pijt=',num2str(Pijt),'Uijt',num2str(U_t),'_092220','.txt'],'wt'); %creat & write Uijt to txt
%      for i=1:i_up
%          for j=1:j_up
%              if j == j_up
%                  fprintf(fid1,'%.4f\n',double(U_count(i,j,U_t))); %0516-keep 4 decimals
%              else
%                  fprintf(fid1,'%.4f\t',double(U_count(i,j,U_t))); %0516-keep 4 decimals
%              end     
%          end    
%      end    
%      fclose(fid1);
% end
    
%   %         
%       fprintf('\nX value: ') 
%       value(X)
%   %             
%   %     fprintf('\nV value: ')
%   %     value(V)
%   %     
%   %     fprintf('\nU_t0 value: ')
%   %     value(U_t0)
%   %     
%   %     fprintf('\nX_t0 value: ') 
%   %     value(X_t0)
%   %     

%   %     fprintf('\nV_t0 value: ')
%   %     value(V_t0)
% 
%   % else
%   % disp('Not solved!')
end

% 092220-keep how 
for U_t=1:t_up 
     fid = fopen(['difCdifP',num2str(P_pa),'UvsD',num2str(U_t),'_101620_v12','.txt'],'wt'); %creat & write Uijt to txt
     for i=1:i_up
         for j=1:j_up
%            if U(i,j,U_t) <= 4.*D_ij_SAV(i,j,U_t)
             if j == j_up
                 fprintf(fid,'%.4f\n',double(U(i,j,U_t)./D_ij_SAV(i,j,U_t))); %0516-keep 4 decimals
%                  fprintf(fid,'%.4f\n',double(U(i,j,U_t)./d(i,j,U_t))); %0516-keep 4 decimals
             else
                 fprintf(fid,'%.4f\t',double(U(i,j,U_t)./D_ij_SAV(i,j,U_t))); %0516-keep 4 decimals
%                  fprintf(fid,'%.4f\t',double(U(i,j,U_t)./d(i,j,U_t))); %0516-keep 4 decimals
             end
%            end
         end    
     end    
     fclose(fid);
end


end

% 0520-save total cost z_op(Pijt) in txt

% fclose(fid);

% 0517-following is save DijSAV in txt
%1.5-save D_ij_SAV to txt for plot use and find Pijt upper bound
% for t=1:t_up 
%     fid = fopen(['D_ij_SAVt',num2str(t),'Pijt=1Iter1.txt'],'wt'); %creat & write D_ij_SAV to txt
%     for i=1:i_up
%         for j=1:j_up
%             if j == j_up
%                 fprintf(fid,'%.4f\n',double(D_ij_SAV(i,j,t))); %0517-keep 4 decimals
%             else
%                 fprintf(fid,'%.4f\t',double(D_ij_SAV(i,j,t))); %0517-keep 4 decimals
%             end     
%         end    
%     end    
%     fclose(fid);
% end



% 0517-following is save X in txt, also include time window one
% 1.7-save Xijt
% for X_SaveTxt_t=1:(t_up+M)
%   fid = fopen(['Xijt',num2str(X_SaveTxt_t),'Pijt=1Iter1.txt'],'wt'); %creat & write Xijt to txt
%   for i=1:i_up
%     for j=1:j_up
%       if j==j_up
%         fprintf(fid,'%.4f\n',double(X(i,j,X_SaveTxt_t))); %0517-keep 4 decimals
%       else
%         fprintf(fid,'%.4f\t',double(X(i,j,X_SaveTxt_t))); %0517-keep 4 decimals
%       end     
%     end    
%    end    
%    fclose(fid);
% end

% 0517-following is save V in txt
% 1.8-save Vit
% fid = fopen(['VitAllPijt=1Iter1.txt'],'wt'); %creat & write Vit to txt
% for V_SaveTxt_t=1:t_up 
%     for i=1:i_up
%       if i==i_up
%         fprintf(fid,'%.4f\n',double(V(i,V_SaveTxt_t))); %0517-keep 4 decimals
%       else
%         fprintf(fid,'%.4f\t',double(V(i,V_SaveTxt_t))); %0517-keep 4 decimals
%       end
%     end        
% end
% fclose(fid);

% 0519-temporarily hide Vi_t0
% 1.9-save Vi_t0 in txt
% fid = fopen(['Vi_t0Pijt=1.txt'],'wt'); %creat & write Vit to txt
% for i=1:i_up
%   fprintf(fid,'%.4f\t',double(V_t0(i,1))); %0517-keep 4 decimals
% end
% fclose(fid);

%} 

%0-end recording running time
toc