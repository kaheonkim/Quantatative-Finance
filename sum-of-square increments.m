clear
clc
randn('state', 100)
 
S0=1;
 
mu=0.05;
sigma=0.3;
L=[100,1000];
t=0;
dt=0.5;
 
for k=1:length(L)
    ddt=dt/L(k);
for j=1:10
    S(j,1)=S0;
    ddS(j,1)=S(j,1)*(mu*ddt+sigma*sqrt(ddt)*randn);
    sum(j,1)=ddS(j,1)^2; 
    SSD(k,1)=S0*sigma^2*dt;
    for i=1:L(k)
        S(j, i+1) = S(j, i) * exp( (mu - 0.5 * sigma * sigma) * ddt + sigma * sqrt(ddt) * randn );
        ddS(j,i+1)=S(j,i+1)*(mu*ddt+sigma*sqrt(ddt)*randn);
        sum(j,i+1)=sum(j,i)+ddS(j,i+1)^2;
        SSD(k,i+1)=S0*sigma^2*dt;
    end
    figure(k)
    plot([t:ddt:t+dt],sum(j,:),'b')
    hold on
 
end
plot([t:ddt:t+dt],SSD(k,:),'r')
title(['The sum-of-square increments with L=',num2str(L(k)),' \newline blue line: sum of square increments, red line:S_0\sigma ^2 \Delta t'])
 
end
 
