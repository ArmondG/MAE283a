close, clear all

importfile1('mass_spring_damper');

%% 3.1
N = length(y);
%G(q) ~ (b1q^-1 +b2q^-2) / (1 + a1q^-1 +a2q^-2)
PHI = [u(2:N-1) u(1:N-2) -y(2:N-1) -y(1:N-2)];
Y = y(3:N);
theta = PHI\Y

%% 3.2
ysim = filter([theta(1) theta(2)],[1 theta(3) theta(4)],u);
plot(t,y)
hold on
plot(t,ysim)
xlabel('time (s)'); ylabel('position (m)')
legend('measured output','simulated output')

%% 3.3
T = 1/10; %From 10Hz sampling freq.
sys = tf([theta(1) theta(2)],[1 theta(3) theta(4)], T);
sysc = d2c(sys, 'zoh');
[num den] = tfdata(sysc, 'v');
newnum = num./num(3);
newden = den./num(3);
%Assuming the zero in numerator is close enough to 0
MassEstimate = newden(1)
DampingEstimate = newden(2)
StiffnessEstimate = newden(3)

