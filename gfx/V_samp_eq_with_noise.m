close all
clear all
%clc

par.K=5.67036713*1e-8;
par.K2=1.38064852*1e-23;    % noise 
par.RTs=800*1e3;
par.ti=65*1e-6;             % integration time
par.tf=1/30;                % frame time
Vb=@(t) Vb_fun(t,par.ti);
par.Vb=@(t) Vb(mod(t,par.tf));
par.Ts=300;
par.C=2.5e-10;              %C=1e-8; 
par.Gleg=2.5e-8;            %1e-7;
timeConst=par.C/par.Gleg;
par.alpha=-0.02;
par.R=@(T) par.RTs*exp(par.alpha*(T-par.Ts));
% we do not know these parameters
par.e=.8;    par.A=(17*1e-6)^2;    par.As=par.A;

% this we know
par.Ps=par.As*par.K*par.Ts^4;

% this we do not know
To=par.Ts;
par.Pt=par.As*par.K*(To+11)^4;
sigma2=1e-4;

% Voltage equation parameters
par.V0=3.1;
par.C2=4*1e-12;
par.E=2;

% euler method
T0=par.Ts;
TT(1)=T0;

rng(0)


sigma=0.005;
output_noise_variance = 1;


N1=100;
N2=100;
M=100;
Vsamp(1)=par.V0;
tt(1)=0;
Vout=[];
for j=1:M
    % solve the problem with Vb>0
    t1=linspace(0,par.ti,N1);
    dt1=t1(2)-t1(1);
    Vsum=0;
    for ii=1:length(t1)-1
        TT(end+1)=TT(end)+dt1*F(t1(ii),TT(end),par,dt1, sigma);
        tt(end+1)=tt(end)+dt1;
%        Vsamp(end+1)=(mod(tt(end),par.ti+par.tf)*par.V0/par.RTs-dt1*sum(par.Vb(tt(end-ii:end))./par.R(TT(end-ii:end))));%+par.E;        

        Vsum=Vsum+par.V0/par.RTs-par.Vb(tt(end))./par.R(TT(end));
        Vsamp(end+1)=(dt1/par.C2)*Vsum+par.E;                         
    end
    Vout(end+1)=Vsamp(end-1);
    

    
    
    % solve the problem with Vb=0
    t2=linspace(par.ti,par.tf,N2);
    dt2=t2(2)-t2(1);
    for ii=1:length(t2)-1
        TT(end+1)=TT(end)+dt2*F(t2(ii),TT(end),par,dt2, sigma);
        tt(end+1)=tt(end)+dt2;
        Vsamp(end+1)=par.E;        
    end

    
end

figure(1); plot(tt,TT,'-m'); grid on
figure(2); plot(tt,Vsamp,'-k'); grid on
figure(3); plot(Vout,'-g')


rng(0)



[TT2, Vout2, Vsamp2]=RunBolometer(par, N1,N2,M, sigma);

figure(1); 
hold on, 
plot(tt,TT2,'-.k'); 
hold off

figure(2); 
hold on
plot(tt,Vsamp2,'-.r');
hold off

figure(3); 
hold on
plot(Vout2,'-.r')
hold off






