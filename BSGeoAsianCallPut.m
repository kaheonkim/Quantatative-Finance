clear
clf
clc
 
S0=10;
E=9;
sigma=0.3;
r=0.05;
T=1;
size=4;
L=12;
dt=T/L;
interval=[0:dt:T];
 
[Call, Put] = BSGeoAsianCallPut_fun(r, S0, E, T, sigma, L);
 
for i=1:size
    M(i)=10^(i+1);
    
    for j=1:M(i)
        S(i,j,1)=1;
        for k=1:L
            S(i,j,k+1)=S0*exp((r-0.5*sigma^2)*dt+sigma*sqrt(dt)*rand);
        end
        SA_mean(i,j)=mean(S(i,j,[2:L+1]));
        SG_mean(i,j)=geomean(S(i,j,[2:L+1]));
        SA_std(i,j)=std(S(i,j,[2:L+1]))/sqrt(L);
        SG_std(i,j)=std(S(i,j,[2:L+1]))/sqrt(L);
    end
    for j=1:M(i)
        VA(i,j)=exp(-r*T)*max(SA_mean(i,j)-E,0);
        VG(i,j)=exp(-r*T)*max(SG_mean(i,j)-E,0);
    end
    
    VA_mean(i)=mean(VA(i,:));
    VG_mean(i)=geomean(VA(i,:));
    
    VA_STD(i)=std(VA(i,:))/sqrt(M(i));
    VG_STD(i)=std(VG(i,:))/sqrt(M(i));
 
end
 
figure(1)
plot([1:size], VA_mean , '--rs')
hold on
plot([1:size], VA_mean - 1.96 * VA_STD, '*b')
hold on
plot([1:size], VA_mean + 1.96 * VA_STD, '*m')
hold on
plot([1:size], Call*ones(1,size), 'LineWidth', 1.5)
figure(2)
plot([1:size],VG_mean , '--rs')
hold on
plot([1:size], VG_mean - 1.96 * VG_STD, '*b')
hold on
plot([1:size], VG_mean + 1.96 * VG_STD, '*m')
hold on
plot([1:size], Call*ones(1,size), 'LineWidth', 1.5)
 
 
function [Call, Put] = BSGeoAsianCallPut_fun(Rate, Price, Strike, Time, sigma, Ave_num)
 
sigma_hat = sigma / Ave_num * sqrt( (Ave_num + 1) * (2 * Ave_num + 1) / 6);
 
mu_hat = (Rate - 0.5 * sigma * sigma) * (Ave_num + 1)/(2 * Ave_num) + 0.5 * sigma_hat * sigma_hat;
 
d1 = (log(Price/Strike) + (mu_hat + 0.5 * sigma_hat * sigma_hat)*Time)/(sigma_hat * sqrt(Time));
 
d2 = d1 - sigma_hat * sqrt(Time);
 
Call = exp(-Rate * Time) * (Price * exp(mu_hat * Time) * normcdf(d1) - Strike * normcdf(d2));
 
Put = exp(-Rate * Time) * (-Price * exp(mu_hat * Time) * normcdf(-d1) + Strike * normcdf(-d2));
end


