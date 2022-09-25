clear
clc
A=1;
S0 = 9;
r = 0.05;
Strike = 10;
sigma = 0.30;
mu = 0.03;
T = 1;
 
for i=1:11
    M(i)=2^(i+9);
    dt=T/M(i);
    
    for j=1:M(i)
        Z=randn;
        dt=T/M(i);
        S(i,j)=S0*exp((r-0.5*sigma^2)*dt*j+sigma*sqrt(T-dt*j)*Z);
        d1(i,j) = (log(S(i,j)/Strike) + (r + 0.5 * sigma * sigma)*(T-dt*j))/(sigma*sqrt(T-dt*j));
        d2 (i,j)= d1(i,j) - sigma * sqrt(T-dt*j);
        Call_Price(i,j) = A*exp(-r*(T-dt*j))*(1-normcdf(d2(i,j)));
        Put_Price(i,j) = A*exp(-r*(T-dt*j))-Call_Price(i,j);
        CDelta(i,j) = A*exp(-r*(T-dt*j))*normpdf(d2(i,j))/(sigma*S(i,j)*sqrt(T-dt*j));
        PDelta(i,j) = -CDelta(i,j);
        if S(i,j)>Strike
            C(i,j)=A;
            P(i,j)=0;
        else
            C(i,j)=0;
            P(i,j)=A;
        end
        C(i,j)=exp(-r*(T-dt*j))*C(i,j);
        P(i,j)=exp(-r*(T-dt*j))*P(i,j);
    end
    C_ave(i)=mean(C(i,1:M(i)));
    P_ave(i)=mean(P(i,1:M(i)));
    C_se(i)=std(C(i,1:M(i)))/sqrt(M(i));
    P_se(i)=std(P(i,1:M(i)))/sqrt(M(i));
    
    for j=1:M(i)-1
        CDE(i,j)=C(i,j+1)-C(i,j);
        PDE(i,j)=P(i,j+1)-P(i,j);
    end
    CDE_ave(i)=mean(CDE(i,1:M(i)-1));
    PDE_ave(i)=mean(PDE(i,1:M(i)-1));
    CDE_se(i)=std(CDE(i,1:M(i)-1))/sqrt(M(i)-1);
    PDE_se(i)=std(PDE(i,1:M(i)-1))/sqrt(M(i)-1);
    
    CDelta_ave(i)=mean(CDelta(i,1:M(i)-1));
    PDelta_ave(i)=mean(PDelta(i,1:M(i)-1));
    CDelta_se(i)=std(CDelta(i,1:M(i)-1))/sqrt(M(i)-1);
    PDelta_se(i)=std(PDelta(i,1:M(i)-1))/sqrt(M(i)-1);
end
MM=[10:20]
figure(1)
plot(MM, C_ave, '--rs')
hold on
plot(MM, C_ave - 1.96 * C_se, '*b')
hold on
plot(MM, C_ave + 1.96 * C_se, '*m')
hold on
xlabel('The number of samples')
legend('Call option value mean', 'Call option value upper interval', 'Call option value lower interval','Black-Scholes-Value')
title('Call Opion Value Approximation(Monte Carlo Method)')
 
figure(2)
plot(MM, P_ave, '--rs')
hold on
plot(MM, P_ave - 1.96 * P_se, '*b')
hold on
plot(MM, P_ave + 1.96 * P_se, '*m')
hold on
 
xlabel('The number of samples')
legend('Put option value mean', 'Put option value upper interval', 'Put option value lower interval','Black-Scholes-Value')
title('Put Opion Value Approximation(Monte Carlo Method)')
figure(3)
plot(MM, CDelta_ave, '--rs')
hold on
plot(MM, CDelta_ave - 1.96 * CDelta_se, '*b')
hold on
plot(MM, CDelta_ave + 1.96 * CDelta_se, '*m')
hold on
 
xlabel('The number of samples')
legend('Call option delta value mean', 'Call option delta value upper interval', 'Call option delta value lower interval','Black-Scholes-Value')
title('Call Opion delta Value Approximation(Monte Carlo Method)')
 
 
figure(4)
plot(MM, PDelta_ave, '--rs')
hold on
plot(MM, PDelta_ave - 1.96 * PDelta_se, '*b')
hold on
plot(MM, PDelta_ave + 1.96 * PDelta_se, '*m')
hold on
xlabel('The number of samples')
legend('Put option delta value mean', 'Put option delta value upper interval', 'Put option delta value lower interval','Black-Scholes-Value')
title('Put Opion delta Value Approximation(Monte Carlo Method)')
 
figure(5)
plot(MM, CDE_ave, '--rs')
hold on
plot(MM, CDE_ave - 1.96 * CDE_se, '*b')
hold on
plot(MM, CDE_ave + 1.96 * CDE_se, '*m')
hold on
xlabel('The number of samples')
legend('Call option delta value mean', 'Call option delta value upper interval', 'Call option delta value lower interval','Black-Scholes-Value')
title('Call Opion delta Value Approximation(Variace Reduction)')
 
 
figure(6)
plot(MM, PDE_ave, '--rs')
hold on
plot(MM, PDE_ave - 1.96 * PDE_se, '*b')
hold on
plot(MM, PDE_ave + 1.96 * PDE_se, '*m')
hold on
xlabel('The number of samples')
legend('Put option delta value mean', 'Put option delta value upper interval', 'Put option delta value lower interval','Black-Scholes-Value')
title('Put Opion delta Value Approximation(Variace Reduction)')
 
 
 
 
 
 


