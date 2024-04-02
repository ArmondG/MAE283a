clear all; close all;

[h,w] = freqz([1,1,1,1,1],[5,0,0,0,0]);

figure(1)
plot(w,abs(h))
xlabel('w [0,pi]'); ylabel('|(H(exp(j*w)|')
%%
close all; clear all;
e = wgn(1000,1,0);
v =[];
for i = 1:1000
    if i < 5
        v(i) = 0;
    else
        v(i) = 1/5*(e(i)+e(i-1)+e(i-2)+e(i-3)+e(i-4));
    end
end
t = 1:1000;
figure(2)
plot(t,v)
hold on
plot(t,e,'.')
xlabel('t');

figure(3)
N = 1000;
t = 0:20;
[c, lags] = xcorr(v,e,20);
stem(lags,c)

    