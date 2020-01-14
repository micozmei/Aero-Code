% Longitudinal Linear LTI State Space Analysis
clear all; close all; clc
%
% DC8 Data. L.V. Schmidt pp. 341-345
%
% Basic Geometry
% ---------------
% S ft^2
S=2600.;
% b ft
b=142.3;
% mac ft
cbar=23.;
%
% Flight Conditions
% -----------------
% Mach number and Altitude (ft)
Mach=0.84;
Alt=33000.;
Theta0=0.;
% speed ft/sec
Speed=825.;
% Dynamic pressure (psf)
qD=271.;
%
% Mass and Inertia
% ----------------
% W lbs
Weight=230000.;
g=32.2; % ft/sec^2
Mass=Weight/g; % slugs
% Iij slug-ft^2
Ixx=3770000.;
Iyy=3560000.;
Izz=7130000.;
Ixz=45000.;
%
% Nondimensional drag derivatives
cD=0.0188;
cDa=0.272;
cDadot=0.;
cDq=0.;
cDM=0.1005;
cDdeltaE=0.0;
%
% Nondimensional lift derivatives
cL=0.326;
cLa=6.744;
cLadot=0.;
cLq=0.;
cLM=0.0;
cLdeltaE=0.352;
%
% Nondimensional pitching moment derivatives
cMa=-2.017;
cMM=-0.17;
cMadot=-6.62;
cMq=-14.60;
cMdeltaE=-1.008;
%
% Dimensional Stability Derivatives (L.V. Schmidt, page 117, Table 4.2)
Xu=-qD*S/Mass/Speed*(2.*cD+Mach*cDM);
Xa=qD*S/Mass*(cL-cDa);
Xadot=-qD*S/Mass*(cbar/2./Speed)*cDadot;
Xq=-qD*S/Mass*(cbar/2./Speed)*cDq;
XdeltaE=-qD*S/Mass*cDdeltaE;
%
Zu=-qD*S/Mass/Speed*(2*cL+Mach*cLM);
Za=-qD*S/Mass*(cD+cLa);
Zadot=-qD*S/Mass*(cbar/2./Speed)*cLadot;
Zq=-qD*S/Mass*(cbar/2./Speed)*cLq;
ZdeltaE=-qD*S/Mass*cLdeltaE;
%
Mu=qD*S*cbar/Iyy/Speed*Mach*cMM;
Ma=qD*S*cbar/Iyy*cMa;
Madot=qD*S*cbar/Iyy*(cbar/2./Speed)*cMadot;
Mq=qD*S*cbar/Iyy*(cbar/2./Speed)*cMq;
MdeltaE=qD*S*cbar/Iyy*cMdeltaE;
%
% The [In], [An], and [Bn] matrices - L.V. Schmidt page 116
In=zeros(4,4);
An=zeros(4,4);
% One input: Elevator
Bn=zeros(4,1);
%
In(1,1)=Speed;
In(2,2)=Speed-Zadot;
In(3,2)=-Madot;
In(3,3)=1.;
In(4,4)=1.;
%
An(1,1)=Speed*Xu;
An(1,2)=Xa;
An(1,4)=-g*cos(Theta0);
An(2,1)=Speed*Zu;
An(2,2)=Za;
An(2,3)=Speed+Zq;
An(2,4)=-g*sin(Theta0);
An(3,1)=Speed*Mu;
An(3,2)=Ma;
An(3,3)=Mq;
An(4,3)=1.;
%
% Elevator effects:
Bn(1,1)=XdeltaE;
Bn(2,1)=ZdeltaE;
Bn(4,1)=MdeltaE;
%
% Find the [A] and [B] matrices for {xdot}=[A]{x}+[B]{u}
A=inv(In)*An;
B=inv(In)*Bn;
%
% Print [A] and compare, in the Lat-Dir case,
% to L.V. Schmidt Example 7.4 page 216
%A
%B
%
% Generate outouts, {y}. made of beta, p, and phi - three outputs
% The matrix [C] transforms the state vector {x} into the output vector {y}
% The following [C] matric is the Unit Matrix. Thus the output vector
% contains four outputs:
% {u/U alpha q Theta }
C=[1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 1. 0.; 0. 0. 0. 1.];
%
% No direct effect of inputs on outputs here:
D=zeros(4,1);
%
% The time over which dynamic response is calculated, and the step size for
% graphing.
t=0.0:0.1:1000.;
%
% Impulse response to elevator command (input 1):
% Note:
%
% If you limit the time integration to 20 sec or so you'll see the short
% period motion getting damped quickly, while the phugoid motion is very
% low-damped with a very long period.
%
% If you integrate over a longer time, say 200 seconds, the short period
% would have been long ago damped, while the phugoid motion with its long
% period and low damping would still continue.
%
[Y,X,t]=impulse(A,B,C,D,1,t);
%
% u/U(t)
Y1=Y(:,1);
% alpha(t):
Y2=Y(:,2);
% q(t):
Y3=Y(:,3);
% theta(t)
Y4=Y(:,4);
% Plot
figure(1)
plot(t,Y1,t,Y2,t,Y3,t,Y4)
title('Impulse')
xlabel('time (s)')
legend('u/U', 'alpha (deg)', 'q (deg/s)', 'theta (deg)', 'Location', 'Best')

eig(A)
figure(2)
hold on
[Numerator, Denominator] = ss2tf(A,B,C,D,1);
H = tf(Numerator(1,:),Denominator);
HH = tf(Numerator(2,:),Denominator);
HHH = tf(Numerator(3,:),Denominator);
HHHH = tf(Numerator(4,:),Denominator);
bode(H) 
bode(HH) 
bode(HHH) 
bode(HHHH)
legend('u/U', 'alpha (deg)', 'q (deg/s)', 'theta (deg)', 'Location', 'Best')