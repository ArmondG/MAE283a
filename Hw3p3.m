clear all; close all; clc;

importfile('R333_STEP_data'); %y and Ts
importfile('R333_FRF_data'); %Gspa and f

N = length(y);

%% 3.1
%Build Hankel Matrices
R    = hankel(y(2:N/2),y(N/2:N-2))-kron(y(1:N/2-1),ones(1,N/2-1));
Rbar = hankel(y(3:N/2+1),y(N/2+1:N-1))-kron(y(2:N/2),ones(1,N/2-1));
[U,S,V]=svd(R);

%Plot the first 20 singular values
s = diag(S);
figure(1)
plot(s(1:20),'*');grid
xlabel('singular value [#]')
ylabel('size of singular value')
title('singular values')

%Type in order of the system
n = 6; %if desired replace 6 with: input("Type the order of the system: ")

%Linear algebra to determine A, B, C, D
R1=U(:,1:n)*sqrt(S(1:n,1:n));R2=sqrt(S(1:n,1:n))*V(:,1:n)';
R1dagger = inv(sqrt(S(1:n,1:n)))*U(:,1:n)';
R2dagger = V(:,1:n)*inv(sqrt(S(1:n,1:n)));
A = R1dagger*Rbar*R2dagger;
B = R2(:,1);
C = R1(1,:);
D = y(1);

%Find the simulated step response
[num_estimate,den_estimate]=tfdata(ss(A,B,C,D),'v');
ysim=step(tf(num_estimate,den_estimate,1),N-1);

%Plot the simulated and original step response
figure(2)
plot(1:N,y,'r',1:N,ysim,'b');
xlabel('time step')
ylabel('output')
title('Comparison of simulated and experimental step response')
legend('system step response','modeled step response')

%% 3.2

%Get initial theta estimates
theta0=[num_estimate,den_estimate];

%Do LS curve fitting
theta_estimate = lsqcurvefit(@(theta, f) myModel(theta, f), ...
    theta0, f, Gspa);

%Plot and compare plots of magnitude and phase angle
Ghat = freqz(theta_estimate(1:7),theta_estimate(8:end),length(f));
figure(3)
plot(f,abs(Ghat),'g-',f,abs(Gspa),'r');
title(['Comparison of Magnitude Response'])
ylabel('|Magnitude|')
xlabel('Frequency (Hz)');
legend('Ghat','Gspa')

figure(4)
plot(f,angle(Ghat),'g-',f,angle(Gspa),'r');
title(['Comparison of Phase Angle'])
ylabel('Phase Angle [rad]')
xlabel('Frequency (Hz)');
legend('Ghat','Gspa')

%Define Model Gmodel
function Gmodel = myModel(theta, freq)
n = length(theta)/2;
b = theta(1:n);
a = theta(n+1:end);
Gmodel = freqz(b,a,length(freq));
end