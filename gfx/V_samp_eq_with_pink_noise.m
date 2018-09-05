close all
clear all
clc

addpath altmany-export_fig-e4117f8/


par.K=5.67036713*1e-8;
par.K2=1.38064852*1e-23;    % noise 
par.RTs=800*1e3;
par.ti=50*1e-6;
par.tf=1/500;
%Vb=@(t) Vb_fun(t,par.ti);
%par.Vb=@(t) Vb(mod(t,par.tf));
par.Vb=3;

par.Ts=300;
par.C=2.5e-10;    %C=1e-8; 
par.Gleg=2.5e-8;   %1e-7;
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
par.V0=5.1;
par.C2=4*1e-12;
par.E=2;

% euler method
T0=par.Ts;
TT(1)=T0;

N1=200;
N2=200;
M=10000;


% generation of the pink noise
NN=1+floor(par.tf/par.ti)*(N1-1);
par.xi = f_alpha(M*NN, 1, 1);

Vsamp(1)=par.V0;
tt(1)=0;
Vout=[];
idx=1;  % used for selecting the pink-noise-vector index
for j=1:M
    % solve the problem with Vb>0
    t1=linspace(0,par.ti,N1);
    dt1=t1(2)-t1(1);
    Vsum=0;
    for ii=1:length(t1)-1
        idx=idx+1;
        TT(end+1)=TT(end)+dt1*F(t1(ii),TT(end),par,dt1,par.Vb,idx);
        tt(end+1)=tt(end)+dt1;

        Vsum=Vsum+par.V0/par.RTs-par.Vb./par.R(TT(end));
        Vsamp(end+1)=(dt1/par.C2)*Vsum+par.E;                         
    end
    Vout(end+1)=Vsamp(end-1);
    

    
    
    % solve the problem with Vb=0
    t2=linspace(par.ti,par.tf,N2);
    dt2=t2(2)-t2(1);
    for ii=1:length(t2)-1
        TT(end+1)=TT(end)+dt2*F(t2(ii),TT(end),par,dt2,0,idx);
        tt(end+1)=tt(end)+dt2;
        Vsamp(end+1)=par.E;        
    end
    idx=idx+(NN-N1+1);
    
end

%figure(1); %clf; hold on; plot(tt,TT,'-m'); hold on; grid on
%figure(2); plot(tt,Vsamp,'-k'); hold on; grid on

%hold on; plot(Vout(10:end),'-g')
plot_spectrum(Vout,1/par.tf)
ylim([1e-11 1e5])
xlim([1 3e2])
legend('Location','southwest')

x_width=15.5;
y_width=6;
set(gcf,'units','centimeters','position',[0 0 x_width y_width])
export_fig -transparent spectrum_pink_noise.eps

function y=F(t,T,par,dt,Vb,idx)
    sigma=(4*par.K2*par.Ts*par.RTs)/(2*par.ti/10);
    sigma=5*sqrt(sigma)*1e-2;
%    sigma=1e-3;
%    sigma=0;
    %sigma=1e-3/6;
    
    y=1/(par.RTs*exp(par.alpha*(T-par.Ts)))...
        *(Vb^2+2*(dt)^(-1/2)*sigma*par.xi(idx)*Vb+par.xi(idx)^2*sigma^2/dt)...
        +par.e*(par.Pt+par.Ps)-(2*par.A)*par.e*par.K*T^4-par.Gleg*(T-par.Ts);
    y=y/par.C;
end