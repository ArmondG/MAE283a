clear, close all; clc;

%% input all the data into MATLAB

%S&P500 data
SaPtable = readtable('SaP'); % taken from finance.yahoo.com
SaPTime = SaPtable.Date';
SaPOpen = SaPtable.Open';
SaPClose = SaPtable.Close';


%Consumer Price Index (CPI)
CPI = [278.802 281.148	283.716	287.504	289.109 292.296	296.311	296.276 ...
    296.171	296.808	298.012	297.711	296.797 299.170	300.840	301.836 ...
    303.363	304.127	305.109	305.691	307.026 307.789	307.671	307.051 307.051];
%taken from www.bls.gov %Date is every month from Dec 2021 to Dec 2023

%Federal Interest Rates
FEDtable = readtable('FED'); %taken from fred.stlouisfed.org
FEDTime = FEDtable.DATE';
FEDIR = FEDtable.FEDFUNDS';
FEDIR(end+1) = FEDIR(end);

%% Get data all on same scale and size

SaPDays = datenum(SaPTime)-datenum('12-15-2021');
FEDDays = datenum(FEDTime)-datenum('12-0-2021'); %Stays the same for whole month including this one
FEDDays(end+1) = datenum('12-0-2023')-datenum('12-0-2021');
CPIDays = FEDDays; %Both on a monthly basis
for i = 1:SaPDays(end)
    Days(i) = i;
    if (Days(i) ~= SaPDays(i))
        SaPClose = [SaPClose(1:i-1) SaPClose(i-1) SaPClose(i:end)];
        SaPOpen = [SaPOpen(1:i-1) SaPOpen(i-1) SaPOpen(i:end)];
        SaPDays = [SaPDays(1:i-1) Days(i-1) SaPDays(i:end)];
    end
    if (Days(i) ~= FEDDays(i))
        FEDIR = [FEDIR(1:i-1) FEDIR(i-1) FEDIR(i:end)];
        CPI = [CPI(1:i-1) CPI(i-1) CPI(i:end)];
        FEDDays = [FEDDays(1:i-1) Days(i-1) FEDDays(i:end)];
        CPIDays = FEDDays;
    end
end

%% Get Auto Correlations of Inputs

N = length(Days);
Ru_SaPOpen = xcorr(SaPOpen, 'biased');
Ru_IR = xcorr(FEDIR, 'biased');
Ru_CPI = xcorr(CPI, 'biased');
tau = -N+1:N-1;

%plot
figure(1)
subplot(3,1,1)
plot(tau,Ru_SaPOpen);
xlabel('tau');ylabel('Ru');
title('S&P 500 Open Price');
subplot(3,1,2)
plot(tau,Ru_IR);
xlabel('tau');ylabel('Ru');
axis([-800, 800, 0, 18]);
title('Interest Rates');
subplot(3,1,3)
plot(tau,Ru_CPI);
xlabel('tau');ylabel('Ru');
title('CPI');
sgtitle(['Auto covariance for N = ' num2str(N) ]);

%% Get Spectral Analysis of Inputs
index=[N:2*N-1 2*N-1 1:N-1];
phiu_SaPOpen = fft(Ru_SaPOpen(index));
phiu_IR = fft(Ru_IR(index));
phiu_CPI = fft(Ru_CPI(index));

% creation of a frequency vector
Deltaf=2*pi/N;
w=[0:Deltaf:pi];

%plot
figure(2)
subplot(3,1,1)
plot(log(w),log(real(phiu_SaPOpen(1:2:N+1))));
xlabel('log(w) [rad/s]');ylabel('log(real{phiu(w)})');
title('S&P 500 Open Price');
subplot(3,1,2)
plot(log(w),log(real(phiu_IR(1:2:N+1))));
xlabel('log(w) [rad/s]');ylabel('log(real{phiu(w)})');
title('Interest Rates');
subplot(3,1,3)
plot(log(w),log(real(phiu_CPI(1:2:N+1))));
xlabel('log(w) [rad/s]');ylabel('log(real{phiu(w)})');
title('CPI');
sgtitle('Estimated Input Spectrum from Auto Covariance');

%Based on the figures above, the frequency range in which the system
%can be reliably estimated appears to be only very low frequencies in the
%ranges from 10^-5 to 10^-2 rad/s

%% SPA estimate of G(e^jw)
G = spa([Days',FEDIR',CPI',SaPOpen'],SaPClose',[],w);
figure(3)
h = bodeplot(G);

%Got the spa estimate from MATLAB directly by plugging in 3 inputs:
%Intersts rates, CPI, and the price open being u1, u2, and u3 respectively.
%Based on the above plots it I believe it will be useful to downsample as
%my system becomes noisier and blows up with higher frequencies.
%By downsampling the model can better estimate the lower frequencies.

%% Downsampling
% Getting a value once every week instead of every day
Index = 1:7:N+1;
Weeks = Days(Index);
SaPOpen = SaPOpen(Index);
SaPClose = SaPClose(Index);
FEDIR = FEDIR(Index);
CPI = CPI(Index);
N = length(Index);

%% Estimate FIR Model

%n < N/10
n = floor(N/10 -1)
i = 3; %number of inputs
u = [CPI' FEDIR' SaPOpen'];
Y = SaPClose(10:N);
for j = 1:3
    PHI = [u(10:N,j) u(9:N-1,j) u(8:N-2,j) u(7:N-3,j) u(6:N-4,j)...
        u(5:N-5,j) u(4:N-6,j) u(3:N-7,j) u(2:N-8,j) u(1:N-9,j)];
    theta(:,j) = PHI\Y';
end

theta

%% Plotting 4 and 5

%Get FIR tf for each input
H = [tf(flip(theta(:,1))',1) tf(flip(theta(:,2))',1) tf(flip(theta(:,2))',1)];

%plot coefficients
k = 1:length(theta);
figure(4)
subplot(3,1,1)
plot(k,theta(:,1),'.r','MarkerSize', 15)
xlabel('k'); ylabel('theta(k)');
title('CPI');
subplot(3,1,2)
plot(k,theta(:,2),'ob','MarkerSize', 10)
xlabel('k'); ylabel('theta(k)');
title('Interest Rates');
subplot(3,1,3)
plot(k,theta(:,3),'k*','MarkerSize', 10)
xlabel('k'); ylabel('theta(k)');
title('S&P 500 Open');
sgtitle('FIR impulse response parameters');


%An important note is with the input of the Interest rate, as seen in
%Figure 3, the magnitude of the G(w) is zero for all frequncies and as such
%doesn't not seem to effect the output system. This is further supported by
%looking at the estimated impulse response coefficients which are on a much
%latrger scale compared to the other two inputs. These large values being
%estimated are the result of compansating for the lack of response by the
%input.

%Get correlation functions
figure(5)
TH1 = arx([SaPClose', CPI'], [0 10 0]);
TH2 = arx([SaPClose', FEDIR'], [0 10 0]);
TH3 = arx([SaPClose', SaPOpen'], [0 10 0]);
%plot
subplot(3,1,1)
resid([SaPClose', CPI'],TH1)
title('CPI','FontSize',8);
subplot(3,1,2)
resid([SaPClose', FEDIR'],TH2)
title('Interest Rates','FontSize',8);
subplot(3,1,3)
resid([SaPClose', SaPOpen'],TH3)
title('S&P 500 Open','FontSize',8)
sgtitle('Estimated Input Spectrum from Auto Covariance', 'FontSize',9);

%% Realization Algorithm
%general I/O data with the generalized realization method

%build Hankel Matrices
y = SaPClose';
N = length(y)-1;
Hy = hankel(y(2:N/2),y(N/2:N-2));
Hybar = hankel(y(3:N/2+1),y(N/2+1:N-1));
Hu1 = hankel(u(2:N/2-1,1), u(N/2+1:N-1,1));
Huo1 = eye(length(Hu1)) - Hu1'*pinv(Hu1*Hu1')*Hu1;
Hu2 = hankel(u(2:N/2-1,2), u(N/2+1:N-1,2));
Huo2 = eye(length(Hu2)) - Hu2'*pinv(Hu2*Hu2')*Hu2;
Hu3 = hankel(u(2:N/2-1,3), u(N/2+1:N-1,3));
Huo3 = eye(length(Hu3)) - Hu3'*pinv(Hu3*Hu3')*Hu3;
H1 = Hy'*Huo1;
H2 = Hy'*Huo2;
H3 = Hy'*Huo3;
Q = [H1 ; H2 ; H3];
[U S V] = svd(Q);

%Plot the singular values
figure(6)
s = diag(S);
plot(s,'*');grid
xlabel('singular value [#]')
ylabel('size of singular value')
title('singular values')

%Based on the graph the system order is 4
%Linear algebra to determine A, B, and C
n = 4;
Q11 = U(1:51,1:n)*sqrt(S(1:n,1:n));
Q12 = sqrt(S(1:n,1:n))*V(:,1:n)';
Q11dagger = inv(sqrt(S(1:n,1:n)))*U(1:51,1:n)';
Q12dagger = V(:,1:n)*inv(sqrt(S(1:n,1:n)));
Q21 = U(51:102,1:n)*sqrt(S(1:n,1:n));
Q22 = sqrt(S(1:n,1:n))*V(:,1:n)';
Q21dagger = inv(sqrt(S(1:n,1:n)))*U(52:102,1:n)';
Q22dagger = V(:,1:n)*inv(sqrt(S(1:n,1:n)));
Q31 = U(103:end,1:n)*sqrt(S(1:n,1:n));
Q32 = sqrt(S(1:n,1:n))*V(:,1:n)';
Q31dagger = inv(sqrt(S(1:n,1:n)))*U(103:end,1:n)';
Q32dagger = V(:,1:n)*inv(sqrt(S(1:n,1:n)));
Qbar1 = Hybar*Huo1; Qbar2 = Hybar*Huo2; Qbar3 = Hybar*Huo2;
Qbar = [Qbar1; Qbar2; Qbar3];
C = [Q11(1,:); Q21(1,:); Q31(1,:)];
A = [Q11dagger'; Q21dagger'; Q31dagger']'*Qbar*[Q12dagger Q22dagger Q32dagger];
B = [Q12(:,1) Q22(:,1) Q32(:,1)];


%Least Squares to solve D
%Build x matrix
x = zeros(105,n,3); %3 is number of inputs
for j = 1:3
    for i = 2:length(x)
        b = A*x(i-1,:)'+C(j,:)'*u(i-1,j);
        x(i,:,j)= b';
    end
end
for j = 1:3
    fun =@(D) y' - [B(:,j)' D]*[x(:,:,j)';u(:,j)'];
    X(j) = lsqnonlin(fun, ones(1,1));
end
D = [X(1) X(2) X(3)];

%% Get impulse response coefficients and plot

for i = 2:n
    g(1,1) = D(1);
    g(i,1) = C(1,:)*A(:,1:4).^(i-1)*B(:,1);
    gss1 = ss(A(:,1:4),B(:,1),C(1,:),D(1));
end
for i = 2:n
    g(1,2) = D(2);
    g(i,2) = C(2,:)*A(:,5:8).^(i-1)*B(:,2);
    gss2 = ss(A(:,5:8),B(:,2),C(2,:),D(2));
end
for i = 2:n
    g(1,3) = D(3);
    g(i,3) = C(3,:)*A(:,9:12).^(i-1)*B(:,3);
    gss3 = ss(A(:,9:12),B(:,3),C(3,:),D(3));
end
figure(7)
subplot(3,1,1)
h1 = impulseplot(gss1);
title('CPI','FontSize',8)
setoptions(h1,'TimeUnits','weeks');
subplot(3,1,2)
h2 = impulseplot(gss2);
title('Interest Rates','FontSize',8)
setoptions(h2,'TimeUnits','weeks');
subplot(3,1,3)
h3 = impulseplot(gss3);
title('S&P 500 Open','FontSize',8)
setoptions(h3,'TimeUnits','weeks');
sgtitle('Impulse Response', 'FontSize',9);

%The accuracy of the FIR model from question 3 is shown to be not very
%accurate. Granted this is a system that is not notoriously difficult to
%model. Nonetheless, the impulse response coefficients found from the
%general realization method are quite different from that from FIR. It can
%be seen that all impulse responses blow up which is expected from the
%system as stock market generally tends upward

%% Residual Plot

figure(8)
subplot(3,1,1)
resid([SaPClose', CPI'],gss1)
title('CPI','FontSize',8)
subplot(3,1,2)
resid([SaPClose', FEDIR'],gss2)
title('Interest Rates','FontSize',8)
subplot(3,1,3)
resid([SaPClose', SaPOpen'],gss3)
title('S&P 500 Open','FontSize',8);
sgtitle('Residuals', 'FontSize',9);

%% Choose Model Structure
% Choosing Box Jenkins because I assume noise is not white and it is easier
% todo multi-input single output estimation of system
nb = [2 2 2];
nc = 1;
nd = 1;
nf = [1 1 1]+1;
nk = [0 0 0];
Ts = 1;
data = iddata(y, [u(:,1), u(:,2), u(:,3)], Ts, 'InputName', ...
    {'CPI', 'Interest Rates', 'S&P 500 Open'}, 'OutputName', 'y');
sys = bj(data, [nb nc nd nf nk]);
sys = chgTimeUnit(sys,'weeks');

%% Plotting the Reu
figure(9)
resid(data,sys)
%Looking at the residuals of the Reu while the values are not exactly zero
%they are still fairly close and are well within the confidence interval,
%likewise the Re appears to have a delta like function which means that the
%error is fairly close to white noise.
%% Plotting Comparison of Y, Ysim, Yprediction
figure(10)
predict(sys, data, 1)
hdat = findobj(gca,'Type','line');
xdat=get(hdat,'Xdata');
ydat=get(hdat,'Ydata');
compare(data, sys);
hold on
plot(xdat{2},ydat{2},'r', 'DisplayName', '1-step Prediction');
ylabel('S&P 500 Stock Price'); xlabel('Time (weeks)');
%Overall the simulation and prediction are quite nice considering the
%dynamics of the system. The chosen order of the system was 2 after some
%iteration this yielded the best simulation and prediction estimates.
%The quality of the simulation and prediction are almost the same, with both 
%having a bit of overshoot over the actual system.







