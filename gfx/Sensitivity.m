function [sens, NETD, noise_std]=Sensitivity(T, sigma, frame_time, integration_time, Frames)

%T=300
%sigma=0.001 
%frame_time=1/30 
%integration_time=65*1e-6



par.K=5.67036713*1e-8;
par.K2=1.38064852*1e-23;    % noise 
par.RTs=800*1e3;

%par.ti=65*1e-6;             % integration time
%par.tf=1/30;                % frame time

par.ti=integration_time;
par.tf=frame_time;

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
par.Pt=par.As*par.K*T^4;

% Voltage equation parameters
par.V0=3.1;
par.C2=4*1e-12;
par.E=2;

N1=100;
N2=100;
M=Frames;

[TT, tt, Vout]=RunBolometer(par, N1,N2,M, sigma);

noise_std=std(Vout(6:end));

M=10;

[TT, tt, Vout2]=RunBolometer(par, N1,N2,M, 0);

par.Pt=par.As*par.K*(T+1)^4;

[TT2, tt2, Vout3]=RunBolometer(par, N1,N2,M, 0);



%figure
%plot(Vout2)
%hold on
%plot(Vout3)
%hold off
%title(['T=' num2str(T) ', frame time=' num2str(frame_time) ', integration time= ' num2str(integration_time)])



sens=abs(Vout3(end)-Vout2(end));


NETD=noise_std/sens;






