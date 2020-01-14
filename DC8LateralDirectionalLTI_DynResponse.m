% Lateral - Directional Linear LTI State Space Analysis
clear all; close all; clc
%
% DC-8 Data. L.V. Schmidt pp. 341-345
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
% Nondimensional Y derivatives
cYbeta=-.7277;
cYp=0.;
cYr=0.;
cYdeltaR=0.1865;
cYdeltaA=0.0;
%
% Nondimensional rolling moment derivatives
clbeta=-0.1673;
clp=-0.516;
clr=0.147;
cldeltaR=0.0211;
cldeltaA=0.0797;
%
% Nondimensional yawing moment derivatives
cnbeta=0.1547;
cnp=-0.011;
cnr=-0.190;
cndeltaR=-0.0834;
cndeltaA=0.0037;
%
% Dimensional Stability Derivatives (L.V. Schmidt, page 117, Table 4.2)
Yb=qD*S/Mass*cYbeta;
Yp=qD*S/Mass*b/2./Speed*cYp;
Yr=qD*S/Mass*b/2./Speed*cYr;
YdeltaR=qD*S/Mass*cYdeltaR;
YdeltaA=qD*S/Mass*cYdeltaA;
%
Lb=qD*S*b/Ixx*clbeta;
Lp=qD*S*b/Ixx*b/2./Speed*clp;
Lr=qD*S*b/Ixx*b/2./Speed*clr;
LdeltaR=qD*S*b/Ixx*cldeltaR;
LdeltaA=qD*S*b/Ixx*cldeltaA;
%
Nb=qD*S*b/Izz*cnbeta;
Np=qD*S*b/Izz*b/2./Speed*cnp;
Nr=qD*S*b/Izz*b/2./Speed*cnr;
NdeltaR=qD*S*b/Izz*cndeltaR;
NdeltaA=qD*S*b/Izz*cndeltaA;
%
% The [In], [An], and [Bn] matrices - L.V. Schmidt page 116
In=zeros(4,4);
An=zeros(4,4);
% Two inputs: Aileron & Rudder
Bn=zeros(4,2);
%
In(1,1)=Speed;
In(2,2)=1.;
In(2,4)=-Ixz/Ixx;
In(3,3)=1.;
In(4,2)=-Ixz/Izz;
In(4,4)=1.;
%
An(1,1)=Yb;
An(1,2)=Yp;
An(1,3)=g*cos(Theta0);
An(1,4)=Yr-Speed;
An(2,1)=Lb;
An(2,2)=Lp;
An(2,4)=Lr;
An(3,2)=1.;
An(4,1)=Nb;
An(4,2)=Np;
An(4,4)=Nr;
%
% Aileron effects:
Bn(1,1)=YdeltaA;
Bn(2,1)=LdeltaA;
Bn(4,1)=NdeltaA;
%
% Rudder Effects:
Bn(1,2)=YdeltaR;
Bn(2,2)=LdeltaR;
Bn(4,2)=NdeltaR;
%
% Find the [A] and [B] matrices for {xdot}=[A]{x}+[B]{u}
A=inv(In)*An;
B=inv(In)*Bn;
%
% Print [A] and compare to L.V. Schmidt Example 7.4 page 216
%A
%B
%
% Generate outouts, {y}. made of beta, p, and phi - three outputs
% The matrix [C] transforms the state vector {x} into the output vector {y}
C=[1. 0. 0. 0. ; 0. 1. 0. 0. ; 0. 0. 1. 0.];
%
% No direct effect of inputs on outputs here:
D=zeros(3,2);
%
% The time over which dynamic response is calculated, and the step size for
% graphing.
t=0.0:0.1:20.;
%
% Impulse response to rudder command (input 2):
[Y,X,t]=impulse(A,B,C,D,2,t);
%
% beta(t)
Y1=Y(:,1);
% p(t):
Y2=Y(:,2);
% phi(t):
Y3=Y(:,3);
% Plot and comapre to L.V. Schmidt Fig. 7.4 page 225
plot(t,Y1,t,Y2,t,Y3)
title('Impulse')
xlabel('time(s)')
legend('Beta (deg)', 'p (deg/s)', 'phi (deg)', 'Location', 'Best')

eig(A)
figure(2)
hold on
[Numerator, Denominator] = ss2tf(A,B,C,D,1);
H = tf(Numerator(1,:),Denominator);
HH = tf(Numerator(2,:),Denominator);
HHH = tf(Numerator(3,:),Denominator);
bode(H) 
bode(HH) 
bode(HHH) 
legend('Beta (deg)', 'p (deg/s)', 'phi (deg)', 'Location', 'Best')