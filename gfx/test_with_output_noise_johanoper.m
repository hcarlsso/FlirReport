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


%sigma=0.000;
%output_noise_variance = 1;


N1=100;
N2=100;
M=100;



par.tf=1/50;                % frame time
par.ti=65*1e-6;             % integration time


Vout={}
TT={}
tt={}

clf


Sigmas = [7*1e-6 logspace(-5,-2,4)]   % Selected noise levels 
T_material= 300 % 50:50:500         % Selected Temperatures


for k=1:length(T_material)   % Iterate over different observed material temperature


%-20+20*k              % Vary Temperature of observed material


%[TT{k}, tt{k}, Vout{k}]=RunBolometer(par, N1,N2,M,0);

%figure(1); 
%hold on, 
%plot(tt{k},TT{k}); 
%hold off

%figure(2); 
%hold on
%plot((1:M),Vout{k})
%hold off
  
% Känslighet
par.Pt=par.As*par.K*(T_material(k)+1)^4;
[TTplus1{k}, ttplus1{k}, Voutplus1{k}]=RunBolometer(par, N1,N2,M,0);
 
par.Pt=par.As*par.K*T_material(k)^4;
[TT{k}, tt{k}, Vout{k}]=RunBolometer(par, N1,N2,M,0);

Denominator= (mean(Voutplus1{k}(100:end))-mean(Vout{k}(100:end)));
    for s=1:length(Sigmas)     % Iterate over different noise levels
       
        [TTnoise{s}, ttnoise{s}, Voutnoise{s}]= ...
            RunBolometer(par, N1,N2,M,Sigmas(s));
       
         NETD(k,s) = std(Voutnoise{s}(100:end))/Denominator

    end
end

   
figure(8)
clf
plot(T_material,log(-NETD))
legend(num2str(Sigmas'))
xlabel('Temperature of observed material (Kelvin)')   
ylabel('log(NETD)')
%legend('50', '100', 'etc')
legend(num2str(T_material'))
