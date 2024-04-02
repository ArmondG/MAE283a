clear all; close all

%% 2.3
N = 1000
a1 = 0.9; 
var = 1;
e = sqrt(3)*(2*rand(N,1)-1);
wt = filter([1 -a1],1,e);
[p, w] = periodogram(wt);
[h, f] = freqz([1 -a1],1,length(w));
phiw = h.^2*var;
plot(log10(w),log10(phiw))
hold on
plot(log10(w),log10(p))
xlabel('log(w)'); ylabel('log(phi(w))')
legend('Expression Based Phi','Periodogram Estimate Phi')
