function [TT, tt,  Vout, Vsamp]=RunBolometer2(par, N1,N2,M, sigma)%, noise_vector)

T0=par.Ts;
TT(1)=T0;

Vsamp(1)=par.V0;
tt(1)=0;
Vout=[];
for j=1:M
    % solve the problem with Vb>0
    t1=linspace(0,par.ti,N1);
    dt1=t1(2)-t1(1);
    Vsum=0;
    for ii=1:length(t1)-1
        TT(end+1)=TT(end)+dt1*F(t1(ii),TT(end),par,dt1, 0); %noisevector
        tt(end+1)=tt(end)+dt1;
%        Vsamp(end+1)=(mod(tt(end),par.ti+par.tf)*par.V0/par.RTs-dt1*sum(par.Vb(tt(end-ii:end))./par.R(TT(end-ii:end))));%+par.E;        

        Vsum=Vsum+par.V0/par.RTs-par.Vb(tt(end))./par.R(TT(end));
        Vsamp(end+1)=(dt1/par.C2)*Vsum+par.E+par.ti*sigma*randn;                         
    end
    Vout(end+1)=Vsamp(end-1);
    

    
    
    % solve the problem with Vb=0
    t2=linspace(par.ti,par.tf,N2);
    dt2=t2(2)-t2(1);
    for ii=1:length(t2)-1
        TT(end+1)=TT(end)+dt2*F(t2(ii),TT(end),par,dt2, 0);
        tt(end+1)=tt(end)+dt2;
        Vsamp(end+1)=par.E;        
    end

    
end





