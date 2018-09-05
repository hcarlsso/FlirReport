

T=300
sigma=0.0001 

frame_time=1/50 
integration_time=65*1e-6
Frames=500


sens=[];
NETD=[];
xgrid=[];
ygrid=[];


sigma_grid=logspace(-2,4, 10)


for k=1:length(sigma_grid)

sigma=sigma_grid(k);

[sens(k), NETD(k)]=Sensitivity2(T, sigma, frame_time, integration_time, Frames)

end

figure(1)

loglog(sigma_grid, NETD)
title('NETD, T=300')
xlabel('sigma')
grid on

% select so that we get 20 mK

sigma=20

Sensitivity2(T, sigma, frame_time, integration_time, Frames)

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

figure(11)

plot(int_times, sens)
title('Sensitivity (Response)')
legend(num2str(frame_times(1)),num2str(frame_times(2)), num2str(frame_times(3)))
grid on

figure(12)
semilogx(int_times, NETD)
title('NETD')
legend(num2str(frame_times(1)),num2str(frame_times(2)), num2str(frame_times(3)))
grid on

figure(13)
plot(int_times, noise_std)
title('noise std')
legend(num2str(frame_times(1)),num2str(frame_times(2)), num2str(frame_times(3)))
grid on



%%

T=300 ;
sigma=0.00001;
Frames=500;

Sensitivity(T, sigma, frame_time, integration_time, Frames)


%1) plot NETD vs sigma, v?lj NETD=0.01-0.03

%2) med ovanst[ende sigma, plotta k'nslighet, NETD som function av frame rate och integretingstid.






