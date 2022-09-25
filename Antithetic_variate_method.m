clear
clc
clf
randn('state', 100) 

Rate = 0.05;
Price = 5;
E = 5.5;
Time = 1;
sigma = 0.5;
Ave_num = 12;
dt = Time/Ave_num;
itr = 7;
A=1;

for i=1:itr
   M(i)=2^(10+i);
   for j=1:M(i)
       S1(1) = Price;
       S2(1) = Price;
       S_anti(1) = 0.5*(S1(1)+S2(1));
       for l=1:Ave_num
           Z = randn;
           S1(l+1) =  S1(l) * exp((Rate - 0.5 * sigma * sigma ) * dt + sigma * sqrt(dt) * Z);
           S2(l+1) =  S2(l) * exp((Rate - 0.5 * sigma * sigma ) * dt + sigma * sqrt(dt) * (-Z));
       end
       
       S1_bar(j) = sum(S1(2:end))/Ave_num;
       
       if S1_bar(j)>E
           V1(j) = exp(-Rate * Time) * A;
       else
           V1(j)=0;
       end
       
       S2_bar(j) = sum(S2(2:end))/Ave_num;
       
       if S2_bar(j)>E
           V2(j) = exp(-Rate * Time) * A;
       else
           V2(j)=0;
       end
       
       
       V_anti(j) = 0.5*(V1(j)+V2(j));
    
   end

    aM(i) = mean(V1);
    bM(i) = std(V1);
    conf_min(i) = aM(i) - 1.96*bM(i)/sqrt(M(i));
    conf_max(i) = aM(i) + 1.96*bM(i)/sqrt(M(i));
    
    %%% Antithetic Estimator
    aM_anti(i) = mean(V_anti);
    bM_anti(i) = std(V_anti);
    conf_anti_min(i) = aM_anti(i) - 1.96 * bM_anti(i)/sqrt(M(i));
    conf_anti_max(i) = aM_anti(i) + 1.96 * bM_anti(i)/sqrt(M(i));   
    
    Ratio(i) = bM(i)/bM_anti(i);
    ITR(i) = 10+i;
end

ITR=ITR';
conf_min=conf_min';
conf_max=conf_max';
conf_anti_min=conf_anti_min';
conf_anti_max=conf_anti_max';
Ratio=Ratio';

table(ITR,conf_min,conf_max,conf_anti_min,conf_anti_max, Ratio)