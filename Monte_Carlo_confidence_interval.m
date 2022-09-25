
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
       S(1) = Price;
       for l=1:Ave_num
           Z = randn;
           S(l+1) =  S(l) * exp((Rate - 0.5 * sigma * sigma ) * dt + sigma * sqrt(dt) * Z);
           
       end
      
       S_bar(j) = sum(S(2:end))/Ave_num;
       if S_bar(j)>E
           V_ar(j) = exp(-Rate * Time) * A;
       else
           V_ar(j)=0;
       end

       
   end

   V_ar_ave(i) = mean(V_ar);
   V_ar_se(i) = std(V_ar)/sqrt(M(i));
   V_LB(i)=V_ar_ave(i)-1.96*V_ar_se(i);
   V_UB(i)=V_ar_ave(i)+1.96*V_ar_se(i);
   ITR(i)=i+10; 
end

V_LB=V_LB';
V_UB=V_UB';
ITR=ITR';

table(ITR,V_LB,V_UB)

