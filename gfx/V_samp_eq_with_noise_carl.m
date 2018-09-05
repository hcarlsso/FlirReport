close all
clear all
clc

addpath altmany-export_fig-e4117f8/


par.K=5.67036713*1e-8;
par.K2=1.38064852*1e-23;    % noise 
par.RTs=800*1e3;
par.ti=65*1e-6;
par.tf=1/30;

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
par.Pts=linspace(par.Pt/2, par.Pt*10,100);
sigma2=1e-4;

% Voltage equation parameters
par.V0=3.1;
par.C2=4*1e-12;
par.E=2;

% euler method
T0=par.Ts;


N1=100;
N2=100;
M=5;



Vb=@(t) Vb_fun_const(t,par.ti); par.Vb=@(t) Vb(mod(t,par.tf));
tt=[]; TT=[]; TT(1)=T0; tt(1)=0; Vout=[]; VoutSamps = []; Vsamp=[]; Vsamp(1)=par.V0;
for pt=1:length(par.Pts)
    par.Pt=par.Pts(pt);
    for j=1:M
        t1=linspace(0,par.ti,N1);
        dt1=t1(2)-t1(1);
        Vsum=0;
        for ii=1:length(t1)-1
            TT(end+1)=TT(end)+dt1*F(t1(ii),TT(end),par,dt1); tt(end+1)=tt(end)+dt1; Vsum=Vsum+par.V0/par.RTs-par.Vb(tt(end))./par.R(TT(end)); Vsamp(end+1)=(dt1/par.C2)*Vsum+par.E;                         
        end
        Vout(end+1)=Vsamp(end-1); t2=linspace(par.ti,par.tf,N2); dt2=t2(2)-t2(1);
        for ii=1:length(t2)-1
            TT(end+1)=TT(end)+dt2*F(t2(ii),TT(end),par,dt2); tt(end+1)=tt(end)+dt2; Vsamp(end+1)=par.E;        
        end
    end
    VoutSamps(pt)=Vout(end);
end
figure(1); hold on;  
AA=plot(par.Pts, VoutSamps,'-g');


Vb=@(t) Vb_fun_triang(t,par.ti); par.Vb=@(t) Vb(mod(t,par.tf));
tt=[]; TT=[]; TT(1)=T0; tt(1)=0; Vout=[]; VoutSamps = []; Vsamp=[]; Vsamp(1)=par.V0;
for pt=1:length(par.Pts)
    par.Pt=par.Pts(pt);
    for j=1:M
        t1=linspace(0,par.ti,N1);
        dt1=t1(2)-t1(1);
        Vsum=0;
        for ii=1:length(t1)-1
            TT(end+1)=TT(end)+dt1*F(t1(ii),TT(end),par,dt1); tt(end+1)=tt(end)+dt1; Vsum=Vsum+par.V0/par.RTs-par.Vb(tt(end))./par.R(TT(end)); Vsamp(end+1)=(dt1/par.C2)*Vsum+par.E;                         
        end
        Vout(end+1)=Vsamp(end-1); t2=linspace(par.ti,par.tf,N2); dt2=t2(2)-t2(1);
        for ii=1:length(t2)-1
            TT(end+1)=TT(end)+dt2*F(t2(ii),TT(end),par,dt2); tt(end+1)=tt(end)+dt2; Vsamp(end+1)=par.E;        
        end
    end
    VoutSamps(pt)=Vout(end);
end
figure(1); hold on;  
BB=plot(par.Pts, VoutSamps,'-b');

Vb=@(t) Vb_fun_bell(t,par.ti); par.Vb=@(t) Vb(mod(t,par.tf));
tt=[]; TT=[]; TT(1)=T0; tt(1)=0; Vout=[]; VoutSamps = []; Vsamp=[]; Vsamp(1)=par.V0;
for pt=1:length(par.Pts)
    par.Pt=par.Pts(pt);
    for j=1:M
        t1=linspace(0,par.ti,N1);
        dt1=t1(2)-t1(1);
        Vsum=0;
        for ii=1:length(t1)-1
            TT(end+1)=TT(end)+dt1*F(t1(ii),TT(end),par,dt1); tt(end+1)=tt(end)+dt1; Vsum=Vsum+par.V0/par.RTs-par.Vb(tt(end))./par.R(TT(end)); Vsamp(end+1)=(dt1/par.C2)*Vsum+par.E;                         
        end
        Vout(end+1)=Vsamp(end-1); t2=linspace(par.ti,par.tf,N2); dt2=t2(2)-t2(1);
        for ii=1:length(t2)-1
            TT(end+1)=TT(end)+dt2*F(t2(ii),TT(end),par,dt2); tt(end+1)=tt(end)+dt2; Vsamp(end+1)=par.E;        
        end
    end
    VoutSamps(pt)=Vout(end);
end
figure(1); hold on;  
CC=plot(par.Pts, VoutSamps,'-m');

Vb=@(t) Vb_fun_lin(t,par.ti); par.Vb=@(t) Vb(mod(t,par.tf));
tt=[]; TT=[]; TT(1)=T0; tt(1)=0; Vout=[]; VoutSamps = []; Vsamp=[]; Vsamp(1)=par.V0;
for pt=1:length(par.Pts)
    par.Pt=par.Pts(pt);
    for j=1:M
        t1=linspace(0,par.ti,N1);
        dt1=t1(2)-t1(1);
        Vsum=0;
        for ii=1:length(t1)-1
            TT(end+1)=TT(end)+dt1*F(t1(ii),TT(end),par,dt1); tt(end+1)=tt(end)+dt1; Vsum=Vsum+par.V0/par.RTs-par.Vb(tt(end))./par.R(TT(end)); Vsamp(end+1)=(dt1/par.C2)*Vsum+par.E;                         
        end
        Vout(end+1)=Vsamp(end-1); t2=linspace(par.ti,par.tf,N2); dt2=t2(2)-t2(1);
        for ii=1:length(t2)-1
            TT(end+1)=TT(end)+dt2*F(t2(ii),TT(end),par,dt2); tt(end+1)=tt(end)+dt2; Vsamp(end+1)=par.E;        
        end
    end
    VoutSamps(pt)=Vout(end);
end
figure(1); hold on;  
DD=plot(par.Pts, VoutSamps,'-k');

EE=plot(par.Pts, 20+5e5*VoutSamps(end)*par.Pts,'--r');
legend([AA BB CC DD EE],'constant','triangular','bell-curve','linear','reference linear slope')
legend('Location','southwest')


xlabel("P_t (m^2 K^4)"); ylabel("V_{samp} (V)"); 
x_width=15; y_width=8; 
set(gcf,'units','centimeters','position',[0 0 x_width y_width])
export_fig -transparent out_vs_inrad.eps


function V=Vb_fun_lin(t,ti)
    beta=3^(3/2)/ti;
    if t<ti
        V=beta*t;
    else
        V=0;
    end
end

function V=Vb_fun_bell(t,ti)

      if t<ti
          mu=ti/2; sigma=ti/10;
          V=3*sqrt(ti)*(2*pi*sigma^2)^(1/4)*normpdf(t, mu, sigma);
      else
          V=0;
      end

end

function V=Vb_fun_triang(t,ti)

    beta=1.5988*1e+05;%7.8653*1e+14;
    alpha=ti*beta;
    if t<ti/2
        V=beta*t;
    elseif ti/2<=t && t<ti
        V=alpha-beta*t;
    else
        V=0;
    end

end

function V=Vb_fun_const(t,ti)
    if t<ti
        V=3*ones(1,length(t));
    else
        V=0*ones(1,length(t));
    end


end


function y=F(t,T,par,dt)
    sigma=(4*par.K2*par.Ts*par.RTs);
    sigma=sqrt(sigma);
    sigma=0;
    %sigma=1e-3;
    
    %sigma=1e-3/6;
    
    y=1/(par.RTs*exp(par.alpha*(T-par.Ts)))...
        *(par.Vb(t)^2+2*(dt)^(-1/2)*sigma*randn*par.Vb(t)+sigma^2)...
        +par.e*(par.Pt+par.Ps)-(2*par.A)*par.e*par.K*T^4-par.Gleg*(T-par.Ts);
    y=y/par.C;
end