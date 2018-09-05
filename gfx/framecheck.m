close all
clear all
clc

addpath altmany-export_fig-e4117f8/


SS=[];
TF =[0.002 0.0035 0.0065 0.0110 0.020];
%TF=0.002;
k=1;
for tf = TF
  par.tf=tf;
  bolom
  SS=[SS; Vout];
  tf
%  Vout=SS(k,:); k=k+1;
  
  sig = Vout(10:end); %remove first samples
  sig = sig-mean(sig);
  q=pwelch(sig);
  n = length(q);
  f=(0:n-1)/(n-1)/tf;
  loglog(f,q); hold on
  grid
  drawnow
end


fmin =min(1./((n-1)*TF));
fmax = max(1./TF);
f = linspace(fmin,fmax,100);
loglog(f,1./f.^2); hold on

hold off
h=legend(num2str(TF')); 
%title('Power spectrum with different frame-rates')
xlabel('frequency (Hz)')
ylabel('S(f)')
axis([fmin fmax 1e-5 1])

x_width=15;
y_width=7;
set(gcf,'units','centimeters','position',[0 0 x_width y_width])
export_fig -transparent pspec.eps
