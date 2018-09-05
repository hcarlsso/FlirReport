function plot_spectrum(sig,fs)

%fs: sampling frequency
sig = sig(100:end); %remove first samples
sig_zeromean=sig-mean(sig);
sig_d=detrend(sig);

q=pwelch(sig,fs);
q_zm=pwelch(sig_zeromean,fs);

q_d=pwelch(sig_d,fs);

j=1:length(q);

figure(1)
title('signal power spectrum using pwelsh')
loglog(j,q,j,q_zm,j,q_d,j,1./j,j,1./j.^2)
grid
legend('signal original','signal zero mean removed','signal detrended','1/f reference','1/f^2 reference')

% 
% 
% figure(2)
% subplot(311)
% periodogram(sig)
% title('periodogram original signal')
% subplot(312)
% periodogram(sig_zeromean)
% title('periodogram zero mean removed')
% subplot(313)
% periodogram(sig_d)
% title('periodogram detrended')
% 
% 
