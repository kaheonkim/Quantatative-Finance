clc
clear
clf
randn('state',100)


S0= 1;
mu=0.05;
sigma=0.5;
T=1;
dt=0.01;
L=T/dt;

MG=[50,100,1000];

for k=1:length(MG)
    
    M=MG(k);
    interval=[0:dt:T];
    for i=1:M
        S(i,1)=S0;
        for j=1:L
            Z=randn(1,L);
            S(i,j+1)=S(i,j)*exp((mu-0.5*sigma^2)*dt+sigma*sqrt(dt)*Z(j));
        end
        SF(i)=S(i,L+1);
    end
    
   

    %% Kernel density function
    dx=0.05;
    center=[0:dx:44*dx];
    
    % Lognormal distribution
    LMSF=mean(log(SF));
    LSSF=sqrt(var(log(SF))); 
    
    %Plot
    N=hist(SF,center);
    figure(2*k-1)
    bar(center, N/(M * dx), 'blue')
    hold on
    plot(center, lognpdf(center,LMSF,LSSF) ,'r','LineWidth',1.5)
    legend('Empirical Distribution', 'logNormal')
    title(['Kernel density estimation (M=',num2str(M),')'])
    
    %% QQ plot

    % Normalization
    MSF=mean(SF);
    SSF=sqrt(var(SF));
    NSF=(SF-MSF)/SSF;
    SNSF = sort(NSF);
    
    % Quantile of normal distribution
    for l = 1:M
        p(l) = l/(1+M);
        Theoretical(l) = norminv(p(l), 0 , 1);
    end
    
    % plot
    figure(2*k)
    scatter(SNSF, Theoretical, 'cyan')
    hold on
    LIN=linspace(min(SNSF),max(SNSF));
    plot(LIN, LIN,'r','LineWidth',1.5)
    legend('Pairs of sorted data and quantile points', 'y=x line')
    title(['QQ plot (M=',num2str(M),')'])

end