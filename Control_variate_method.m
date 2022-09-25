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

sigma_hat = sigma / Ave_num * sqrt( (Ave_num + 1) * (2 * Ave_num + 1) / 6);

mu_hat = (Rate - 0.5 * sigma * sigma) * (Ave_num + 1)/(2 * Ave_num) + 0.5 * sigma_hat * sigma_hat;
d1 = (log(Price/E) + (mu_hat + 0.5 * sigma_hat * sigma_hat)*Time)/(sigma_hat * sqrt(Time));
d2 = d1 - sigma_hat * sqrt(Time);
Call=A*exp(-Rate*Time)*normcdf(d2);


for i=1:itr
    M(i)=2^(10+i);
    for j=1:M(i)
        S(1) = Price;
        for l=1:Ave_num
           Z = randn;
           S(l+1) =  S(l) * exp((Rate - 0.5 * sigma * sigma ) * dt + sigma * sqrt(dt) * Z);
           
        end
        S_geomean(j)=geomean(S(2:end));
        
        if S_geomean(j)>E
           V_geo(j) = exp(-Rate * Time) * A;
        else
            V_geo(j)=0;
        end
        
        S_arimean(j)=sum(S(2:end))/Ave_num;
        
        if S_arimean(j)>E
           V_ari(j) = exp(-Rate * Time) * A;
        else
            V_ari(j)=0;
        end
        
    end
    
    V_=V_ari+Call-V_geo;
    
    %%% Standard Estimator
    aM(i) = mean(V_ari);
    bM(i) = std(V_ari);
    conf_min(i) = aM(i) - 1.96*bM(i)/sqrt(M(i));
    conf_max(i) = aM(i) + 1.96*bM(i)/sqrt(M(i));
    
    %%% Control Variate
    aM_cv(i) = mean(V_);
    bM_cv(i) = std(V_);
    conf_cv_min(i) = aM_cv(i) - 1.96 * bM_cv(i)/sqrt(M(i));
    conf_cv_max(i) = aM_cv(i) + 1.96 * bM_cv(i)/sqrt(M(i));   
    
    Ratio(i) = bM(i)/bM_cv(i);
    ITR(i) = 10+i;
    
end
ITR=ITR';
conf_min=conf_min';
conf_max=conf_max';
conf_cv_min=conf_cv_min';
conf_cv_max=conf_cv_max';
Ratio=Ratio';

table(ITR,conf_min,conf_max,conf_cv_min,conf_cv_max, Ratio)
