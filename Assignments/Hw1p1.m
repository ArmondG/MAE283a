clear all; close all;
importfile('noisy_sine_data')
N =1000;
dT =1/1000;
samplingfreq = 1/dT;

freqresolution = samplingfreq/N
maxfreq = 0.5*samplingfreq

dw = 2*pi/(N*dT);
Un = [];
for w = 1:maxfreq
    U = 0;
    for k = 1:N
         U = U + u(k)*exp(j*w*dw*k*dT);
    end
    Un(w) = U;
end
w = 1:freqresolution:maxfreq;
aUn = (abs(Un));
plot(w,aUn)
xlabel('frequency Hz'); ylabel('|Un(wm)|')
