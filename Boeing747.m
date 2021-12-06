clc
syms x
%%at SSL with M = 0.25, 20deg flaps extended
M = 0.25; rho=1.225;%kg/m^3 
R=287;%Jkg/K 
gamma=1.4; 
T=288;%K
g = 9.81; %m/s^2
cl0 = 1.11; cd0=0.102; cla=5.7; cda=0.66; cma=-1.26; cmad=-3.2; cmq = -20.8;
u0 = M*sqrt(gamma*R*T); %m/s
Q = 0.5*rho*(u0^2); %Pa
W = 2508939.33; %N
m = W/g;
Iy = 43798800; %kg*m^2
s = 525; %m^2
c = 9; %m

Xu = -2*cd0*Q*s/(u0*m);
Zu = -2*cl0*Q*s/(u0*m);
Xw = -(cda-cl0)*Q*s/(u0*m);
Zw = -(cla+cd0)*Q*s/(u0*m);
Mw = cma*Q*s*c/(u0*Iy);
Mwd = cmad*Q*s*(c^2)/(2*(u0^2)*Iy);
Mq = cmad*Q*s*(c^2)/(2*u0*Iy);

A = [Xu Xw 0 -g;Zu Zw u0 0;Mwd+(Mwd*Zu) Mw+(Mwd*Zw) Mq+(Mwd*u0) 0;0 0 1 0];
Eigenval = det((x*eye(4))-A);
CE = sym2poly(Eigenval);
TF = tf(1,CE);

%phugoid
Phug = [A(1,1) A(1,4);-A(2,1)/u0 A(2,4)];
S1 = det((x*eye(2))-Phug);
C1 = sym2poly(S1);
r1 = roots(C1); %phugoid roots

wnph = sqrt(C1(3));
zetaph = C1(2)/(2*wnph);

%Short-Period
Za = u0*Zw;
Ma = u0*Mw;
Mad = u0*Mwd;

SP = [Zw 1;Ma+(Mad*(Za/u0)) Mq+Mad];
S2 = det((x*eye(2))-SP);
C2 = sym2poly(S2);
r2 = roots(C2); %SP roots

wnsp = sqrt(C2(3));
zetasp = C2(2)/(2*wnsp);

%plotting the roots
r1r = real(r1);
r1i = imag(r1);
r2r = real(r2);
r2i = imag(r2);

figure('Name','Phugoid and short period roots')
plot(r1r(1),r1i(1),'Marker','*','Color',[1 0 0])
hold on
plot(r1r(2),r1i(2),'Marker','*','Color',[1 0 0])
hold on
plot(r2r(1),r2i(1),'Marker','*','Color',[1 0 0])
hold on
plot(r2r(2),r2i(2),'Marker','*','Color',[1 0 0])
grid on
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
xlabel('Real-axis')
ylabel('Imaginary-axis')
hold off

figure('Name','root locus')
rlocus(TF)
