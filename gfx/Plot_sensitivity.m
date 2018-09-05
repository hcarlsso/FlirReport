close all
clear all
clc

addpath altmany-export_fig-e4117f8/


T=300
sigma=0.0001 

frame_time=1/50 
integration_time=65*1e-6
Frames=500


sens=[];
NETD=[];
xgrid=[];
ygrid=[];


sigma_grid=logspace(-7,-4, 5)


for k=1:length(sigma_grid)

sigma=sigma_grid(k);

[sens(k), NETD(k)]=Sensitivity(T, sigma, frame_time, integration_time, Frames)

end

figure(1)
x_width=15; y_width=5; 
loglog(sigma_grid, NETD*1000, 'LineWidth',2)
title('NETD at T=300')
xlabel('\sigma')
ylabel('NETD (mK)')
grid on
set(gcf,'units','centimeters','position',[0 0 x_width y_width])
export_fig -transparent NETD_as_function_of_sigma.eps


% select so that we get 20 mK

sigma=7*1e-6

Sensitivity(T, sigma, frame_time, integration_time, Frames)

%(4*par.K2*par.Ts*par.RTs)


%%

T=300
sigma=7*1e-6

frame_time=1/30 
integration_time=65*1e-6
Frames=100


sens=[];
NETD=[];
noise_std=[];

frame_times=[1/30, 1/50, 1/100]
int_times=linspace(5, 100, 10)*1e-6

for k=1:length(frame_times)
    for ell=1:length(int_times)
        frame_time=frame_times(k); 
        integration_time=int_times(ell);
        
        [sens(k,ell), NETD(k,ell), noise_std(k, ell)]=Sensitivity(T, sigma, frame_time, integration_time, Frames);

    end
    
end

figure(2)
plot(int_times, sens, 'LineWidth',2)
title('Sensitivity (Response)')
ylabel('Sensitivity (V)')
xlabel('Integration time (s)')
legend('Frame rate 30 Hz', 'Frame rate 50 Hz','Frame rate 100 Hz')
grid on

%x_width=5; y_width=5; 
%set(gcf,'units','centimeters','position',[0 0 x_width y_width])
%export_fig -transparent NETD_as_function_of_sigma.eps

figure(3)
plot(int_times, NETD*1000, 'LineWidth',2)
title('NETD')
ylabel('NETD (mK)')
xlabel('Integration time (s)')
legend('Frame rate 30 Hz', 'Frame rate 50 Hz','Frame rate 100 Hz','Location','northwest')
grid on

x_width=8; y_width=7; 
set(gcf,'units','centimeters','position',[0 0 x_width y_width])
export_fig -transparent NETS_Function_of_Integration_Time.eps
y=[1e-5, 1e-4];

figure(4)
plot(int_times, noise_std, 'LineWidth',2)
title('Noise std')
ylabel('Noise std (V)')
xlabel('Integration time (s)')
legend('Frame rate 30 Hz', 'Frame rate 50 Hz','Frame rate 100 Hz','Location','northwest')
grid on
set(gcf,'units','centimeters','position',[0 0 x_width y_width])
export_fig -transparent STD_Function_of_Integration_Time.eps


%%

T=300 ;
sigma=0.00001;
Frames=500;

Sensitivity(T, sigma, frame_time, integration_time, Frames)


%1) plot NETD vs sigma, v?lj NETD=0.01-0.03

%2) med ovanst[ende sigma, plotta k'nslighet, NETD som function av frame rate och integretingstid.






