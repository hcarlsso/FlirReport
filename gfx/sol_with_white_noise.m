close all
clear all
clc
addpath altmany-export_fig-e4117f8/

par.K=5.67036713*1e-8;
par.K2=1.38064852*1e-23;    % noise 
par.RTs=800*1e3;
par.ti=65*1e-6;
par.tf=1/30;
Vb=@(t) Vb_fun(t,par.ti);
par.Vb=@(t) Vb(mod(t,par.tf));
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
par.Pt=par.As*par.K*(To+10)^4;
sigma2=1e-4;

% Voltage equation parameters
par.V0=3.1;
par.C2=4*1e-12;
par.E=2;

% euler method
T0=par.Ts;
TT(1)=T0;
Vsamp(1)=par.V0;
Vout=[];

N1=10000;
N2=10000;
M=3;

tt(1)=0; tout=[];
for j=1:M
    % solve the problem with Vb>0
    t1=linspace(0,par.ti,N1);
    dt1=t1(2)-t1(1);
    Vsum=0;
    for ii=1:length(t1)-1
        TT(end+1)=TT(end)+dt1*F(t1(ii),TT(end),par,dt1,0);
        tt(end+1)=tt(end)+dt1;
        
        Vsum=Vsum+par.V0/par.RTs-par.Vb(tt(end))./par.R(TT(end));
        Vsamp(end+1)=(dt1/par.C2)*Vsum+par.E; 
    end
    Vout(end+1)=Vsamp(end-1);
    tout(end+1)=tt(end-1); 

    
    % solve the problem with Vb=0
    t2=linspace(par.ti,par.tf,N2);
    dt2=t2(2)-t2(1);
    for ii=1:length(t2)-1
        TT(end+1)=TT(end)+dt2*F(t2(ii),TT(end),par,dt2,0);
        tt(end+1)=tt(end)+dt2;
        Vsamp(end+1)=0*par.E;                
    end

    
end

Vout(end+1)=Vout(end);
tout(end+1)=tt(end-1); 

figure(1)
NF=N1+N2;
plot(tt(NF+1:NF+N1-2),TT(NF+1:NF+N1-2),'-k'); hold on; grid on

figure(2)
plot(tt(NF+N1+1:NF+N1+N2),TT(NF+N1+1:NF+N1+N2),'-k'); hold on; grid on


figure(3)
plot(tt,TT,'-k'); hold on; grid on 

figure(4)
plot(tout,Vout,'-*k'); hold on; grid on


%% solution with noise
sigma_w=(4*par.K2*par.Ts*par.RTs);
sigma_w=5e3*sqrt(sigma_w);
tt=[]; TT=[];




% euler method
T0=par.Ts;
TT(1)=T0;
Vsamp(1)=par.V0;
Vout=[];

N1=10000;
N2=10000;
M=3;
tt(1)=0; tout=[];
for j=1:M
    % solve the problem with Vb>0
    t1=linspace(0,par.ti,N1);
    dt1=t1(2)-t1(1);
    Vsum=0;
    for ii=1:length(t1)-1
        TT(end+1)=TT(end)+dt1*F(t1(ii),TT(end),par,dt1,sigma_w);
        tt(end+1)=tt(end)+dt1;
        
        Vsum=Vsum+par.V0/par.RTs-par.Vb(tt(end))./par.R(TT(end));
        Vsamp(end+1)=(dt1/par.C2)*Vsum+par.E; 
    end
    Vout(end+1)=Vsamp(end-1);
    tout(end+1)=tt(end-1); 

    
    % solve the problem with Vb=0
    t2=linspace(par.ti,par.tf,N2);
    dt2=t2(2)-t2(1);
    for ii=1:length(t2)-1
        TT(end+1)=TT(end)+dt2*F(t2(ii),TT(end),par,dt2,sigma_w);
        tt(end+1)=tt(end)+dt2;
        Vsamp(end+1)=0*par.E;                
    end

    
end

Vout(end+1)=Vout(end);
tout(end+1)=tt(end-1); 

x_width=7; y_width=6; 
figure(1)
NF=N1+N2;
plot(tt(NF+1:NF+N1-2),TT(NF+1:NF+N1-2),'--r'); hold on; grid on
xlabel("time (secs)"); ylabel("Temperature (K)"); 
%title("Integration time")
set(gcf,'units','centimeters','position',[0 0 x_width y_width])
%saveas(gcf,'fig1_int_time','jpeg')
export_fig -transparent sol_int_time_noise.eps

figure(2)
plot(tt(NF+N1+1:NF+N1+N2),TT(NF+N1+1:NF+N1+N2),'--r'); hold on; grid on
xlabel("time (secs)"); ylabel("Temperature (K)"); 
%title("Cooling time")
set(gcf,'units','centimeters','position',[0 0 x_width y_width])
%saveas(gcf,'fig1_cooling','jpeg')
export_fig -transparent sol_cooling_noise.eps


x_width=15.5; 
figure(3)
plot(tt,TT,'--r'); grid on 
xlabel("time (secs)"); ylabel("Temperature (K)"); %title("Solution with several pulses")
set(gcf,'units','centimeters','position',[0 0 x_width y_width])
export_fig -transparent sol_several_pulses_noise.eps


figure(4)
plot(tout,Vout,'--or'); grid on 
xlabel("time (secs)"); ylabel("V_{samp} (V)"); 
%title("Readout with several pulses")
set(gcf,'units','centimeters','position',[0 0 x_width y_width])
%saveas(gcf,'fig1_several_pulses','jpeg')
export_fig -transparent Vout_several_pulses_noise.eps





function V=Vb_fun(t,ti)
    VVb=3;
    if t<ti
        %V=VVb-(VVb/ti)*t;
        V=VVb;
    else
        V=0;
    end
end

        
function y=F(t,T,par,dt,sigma_w)

    %sigma=1e-3;
    
    %sigma=1e-3/6;
    
    y=1/(par.RTs*exp(par.alpha*(T-par.Ts)))...
        *(par.Vb(t)^2+2*(dt)^(-1/2)*sigma_w*randn*par.Vb(t)+sigma_w^2)...
        +par.e*(par.Pt+par.Ps)-(2*par.A)*par.e*par.K*T^4-par.Gleg*(T-par.Ts);
    y=y/par.C;
end