% 737
clear; close all; clc
% l=roll, m=pitch, n=yaw, alpha=AOA, Beta=sideslip, Y=sideforce

% \\\\\\\\\\\\\\\\\\\\ Geometric Inputs \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
a_0_w                = 5.9;        %lift curve slope at zero lift for wing (per radian)
a_0_t                = 5.9;        %lift curve slope for horizontal tail at zero lift (per radian)
alpha_v              = 15;          %Angle of attack of vertical tail in degrees
alpha_w              = 15;        %Angle of attack of wing in degrees
b_a                  = 12.437;       %Aileron span
b_t                  = 47.083;       %Horizontal tail span
b_v                  = 25.7;       %Vertical tail span
b_w                  = 117.417;      %Wing span
BRA                  = 0;      %Body reference area (fuselage volume^2/3)
c_a                  = 3.93;        %Aileron chord
c_t                  = 8.51;         %Mean Aerodynamic Chord of tail
C_D_t                = 0;          %Drag coefficient for the tail (3D)
C_l                  = 0.38;       %Aircraft lift coefficient (2D)
C_l_alpha_t          = 5.9;        %Tail lift curve slope (2D)
C_l_alpha_w          = 5.9;        %Wing lift curve slope (2D)
C_L_alpha_f          = 0.1;        %Fuselage lift curve slope (3D)
C_L_f                = 0;          %Lift coefficient for the fuselage (3D)
C_m_ac               = 0;          %Moment coefficient at A.C. of wing
C_m_ac_t             = 0;          %Moment coefficient at A.C. of horizontal tail
C_Y_p_over_C_L       = 0.15;       %per radian
d                    = 12.333;      %Maximum diameter of fuselage
delta_C_l_Beta_1     = 0.0008;     %Wing position correction
delta_C_l_Beta_2     = -0.00016;   %Wing position correction
dC_d_0_over_dalpha   = 0;          %per radian
Delta_C_l_Beta_tip   = 0;          %Wingtip shape correction factor
e                    = 0.99;       %Oswald efficiency factor
e1_t                 = 0.99;       %Oswald efficiency factor of tail
e1_w                 = 0.99;       %Oswald efficiency factor of wing
eta_t                = 0.85;       %(q_t/q)
eta_v                = 0.55;       %(q_v/q)
Gamma                = 8;        %Dihedral angle in degrees
h_t                  = 0;        %Height of the horizontal tail A.C. above the C.G. 
i_t                  = 0;         %Horizontal tail incidence angle
i_w                  = 0;        %Wing incidence angle
k_t                  = 1.25;       %Wingtip correction factor
K_n                  = 0.1455;     %Empirical factor as function of fitness ratio and C.G. location (per radian)
l_b                  = 110.33;        %Length of fuselage
l_t                  = 0;       %Horizontal distance between C.G. and horizontal tail quarter chord
l_t_prime            = 0;       %Horizontal distance from wing quarter chord to horizontal tail quarter chord
l_v                  = 0;       %Horizontal distance between C.G. and A.C. of vertical tail
lambda               = 0.15;       %Taper ratio (c_t/c_r)
Lambda               = 25;         %Wing sweep angle in degrees
Lambda_t             = 35;       %Tail sweep angle in degrees
r_1                  = 10;       %Fuselage diameter in tail region
rho                  = 0.001496;   %Air density
S_e                  = 77.6;        %Surface area of elevator
S_r                  = 147.8;        %Surface area of rudder
S_s                  = 0;        %Body side area
S_t                  = 425.95;        %Surface area of horizontal tail
S_v                  = 347.9;        %Surface area of vertical tail
S_w                  = 1341.2;       %Surface area of wing
u                    = 421.952;      %Airspeed
w_f                  = 0;      %Maximum width of fuselage or nacelles
W                    = 174200;     %Aircraft weight in lbs
x_a                  = 0;       %Distance parallel to relative wind from A.C. to C.G. of wing
x_ac                 = 0;       %Distance from nose to A.C. of wing
x_prime              = 0;      %Distance from C.G. to wing quarter chord
Y_i                  = 34.04;       %Spanwise distance from centerline to the inboard edge of the control surface 
z_a                  = 0;      %Vertical distance from C.G. to A.C. of wing
z_T                  = 0;       %Vertical distance from C.G. to thrust line
z_v                  = 0;       %Vertical distance from x-axis to A.C. of vertical tail
z_w                  = 0;       %Vertical distance from wing root quarter chord to fuselage centerline
A_nacelle            = 20.295; %Area of nacelle
n_nacelle            = 2;      %Number of nacelles
f_w                  = 0.008*S_w;
f_f                  = 0.063*pi*(d/2)^2+0.12*(n_nacelle*A_nacelle);
f_t                  = 0.0043*S_t;
A_gear               = 0;
f_maingear           = 0;
f_nosegear           = 0; 


% \\\\\\\\\\\\\\\\\\\\ Calculations \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
q                       = 0.5*rho*u*u;                                                    %Dynamic pressure of airplane
c                       = 22.9; %S_w/b_w;                                                   %Mean Aerodynamic Chord of wing
AR                      = 9.44; %b_w/c;                                                          %Aspect ratio of wing
Ae                      = 1.55*(b_v*b_v)/S_v;                                             %Effective aspect ratio
AR_t                    = b_t/c_t;                                                        %Aspect ratio of tail
C_L_w                   = C_l/(1+2/AR);                                                   %Lift coefficient for finite wing
C_L_t                   = C_L_w*(x_a/l_t)*(S_w/S_t)*(1/eta_t);                            %Lift coefficient of tail
%C_L_alpha_t            = 57.3*C_l_alpha_t/(1+(57.3*C_l_alpha_t)/(pi*e1_t*AR_t));         %Tail lift curve slope in degrees
%C_L_alpha_w            = 57.3*C_l_alpha_w/(1+(57.3*C_l_alpha_w)/(pi*e1_w*AR));           %Wing lift curve slope in degrees
C_L_alpha_t             = C_l_alpha_t*(AR_t/(2+sqrt(4+AR_t^2)));
C_L_alpha_w             = C_l_alpha_w*(AR/(2+sqrt(4+AR^2)));
C_D_0                   = (f_w+f_f+f_maingear+f_nosegear+f_t)/S_w;                        %Zero lift or parasitic drag (3D)
epsilon                 = 20*C_L_w*((3*c/l_t_prime)^0.25)*((1/lambda)^0.3)/(AR^0.725);    %Downwash angle
depsilon_over_dalpha    = 20*(C_L_alpha_w/57.3)*((3*c/l_t_prime)^0.25)*((1/lambda)^0.3)/(AR^0.725); %Change in downwash angle due to change in angle of attack in radians
alpha_t                 = alpha_w-i_w+i_t-epsilon;                                        %Angle of attack of tail
alpha_required          = C_L_t/(C_L_alpha_t/57.3);                                       %Angle of attack to achieve lift coefficient
k_f                     = figure9(x_ac,l_b);                                              %Fuselage correction factor
dalpha_t_over_ddelta_e  = figure13(S_e,S_t);                                              %per radian
delta_e                 = (alpha_required-alpha_t)/dalpha_t_over_ddelta_e;                %Elevator deflection
K_i                     = figure23(z_w,d);                                                %Wing-fuselage interference factor
k                       = figure24(b_v,r_1);                                              %Function of ratio of vertical tail span to fuselage diameter in tail region
C_l_Beta_over_Gamma     = figure26(AR,lambda);                                            %per degree
a_v                     = figure28(Ae)*57.3;                                              %Lift curve slope of the vertical tail (per radian)
Delta_C_Y_p_Gamma_over_C_l_p_Gamma0 = figure33(Gamma);                                    %per radian
C_n_p_over_C_L_Lambda0  = figure35(lambda,AR);                                            %per radian
tau2                    = figure45(S_r,S_v);                                              %Function of rudder area to vertical tail area ratio
dsigma1                 = figure36(lambda,z_v,b_w);                                       %Effect of wing sidewash                                                          
dsigma2                 = 9.3*(3/b_w)*((z_v-(z_v*cosd(alpha_w)-l_v*sind(alpha_w)))/b_w)^2;%Effect of wing sidewash
C_l_p_a_0w              = figure34(lambda,AR);                                            %per radian
C_l_p_a_0t              = figure34(lambda,AR_t);                                          %per radian
C_L_delta_A_over_taui   = figure42(lambda,Y_i,0,b_w);                                     %per radian
C_L_delta_A_over_tauo   = figure42(lambda,Y_i,b_a,b_w);                                   %per radian
tau1                    = figure41(c_a,c);                                                %Function of aileron chord to wing chord ratio 
K                       = figure43(Y_i,b_w,AR)+0.05;                                      %Empirical factor as function of semi-span and spanwise distance from center to inboard (for taper ratio = 0.5), added 0.05 as correction factor
C_l_delta_e             = 57.3*figure14(c_a,c_t,AR_t);
%delta_curve = polyfit([0 0.0625 0.125 0.375 0.625 0.875 1],[0.14 0.08 0.045 0.01 0.02 0.04 0.05],3);
%delta_w = polyval(delta_curve, lambda);
%e = 1/(1 + delta_w);


% \\\\\\\\\\\\\\\ Stability Derivatives \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
% Longitudinal
C_L               = W/(q*S_w); %C_L_w+C_L_t*(S_t/S_w)*(eta_t)+C_L_f;
C_D               = C_D_0+((C_L)^2)/(pi*e*AR);
T                 = C_D*q*S_w; %Thrust
C_m               = C_L*(x_a/c)+C_D*(z_a/c)+C_m_ac-C_L_t*(S_t/S_w)*(l_t/c)*(eta_t)+C_D_t*(S_t/S_w)*(h_t/c)*(eta_t)+C_m_ac_t;
C_T               = T/(0.5*rho*S_w*u^2);
C_m_T             = (T*z_T)/(q*S_w*c);
C_L_u             = 0;
C_D_u             = 0;
C_m_u             = 0;
C_L_alpha         = C_L_alpha_w+C_L_alpha_f+(C_L_alpha_t)*(1-depsilon_over_dalpha)*(S_t/S_w)*(eta_t);
C_D_alpha         = dC_d_0_over_dalpha+(2*C_L/(pi*e*AR))*C_L_alpha;
C_m_alpha_w       = ((1+((2*C_L)/(pi*e*AR))*(alpha_w-i_w)+C_D/C_L_alpha)*(x_a/c)+((2*C_L)/(pi*e*AR)-(alpha_w-i_w)-(C_L/C_L_alpha))*(z_a/c))*C_L_alpha;
C_m_alpha_f       = k_f*w_f*w_f*l_b/(S_w*c);
C_m_alpha_t       = (C_L_alpha_t)*(1-depsilon_over_dalpha)*(S_t/S_w)*(l_t/c)*(eta_t);
C_m_alpha         = (C_m_alpha_w+C_m_alpha_f-C_m_alpha_t)*(-0.5); %Correction Factor Added
C_L_alphadot      = 2*(C_L_alpha_t)*depsilon_over_dalpha*(l_t_prime/c)*(S_t/S_w)*(eta_t);
C_D_alphadot      = 0;
C_m_alphadot      = (-2*C_L_alpha_t*depsilon_over_dalpha*(l_t/c)*(l_t_prime/c)*(S_t/S_w)*(eta_t))*(0.84); %Correction Factor Added
C_L_q             = 2*(x_prime/c)*C_L_alpha+2*(l_t/c)*C_L_alpha_t*(S_t/S_w)*(eta_t);
C_D_q             = 0;
C_m_q             = (((-2*x_prime/c^2)*abs(x_prime)*C_L_alpha-((2*(l_t)^2)/c^2)*C_L_alpha_t*(S_t/S_w)*(eta_t)))*(2); %Correction Factor Added
C_L_delta_e       = 1.05*C_l_delta_e*(C_L_alpha_t/(C_l_alpha_t))*(S_t/S_w)*(eta_t);
C_D_delta_e       = 0*(S_t/S_w)*eta_t;
C_m_delta_e       = (-l_t/c)*C_L_delta_e;

% Lateral
C_Y_Beta_t        = -k*(a_v)*(0.724+1.53*(S_v/S_w)+0.4*(z_w/d)+0.009*AR)*(S_v/S_w);
C_Y_Beta_f        = -K_i*C_L_alpha_f*(BRA/S_w);
C_Y_Beta_w        = -0.0001*abs(Gamma)*57.3;
C_Y_Beta          = (C_Y_Beta_t+C_Y_Beta_f+C_Y_Beta_w)*0.5; %Correction Factor
C_l_Beta_w        = (C_l_Beta_over_Gamma)*Gamma*57.3+Delta_C_l_Beta_tip;
C_l_Beta_w_Gamma0 = C_L*(0.05-k_t*(0.29+0.71*lambda)/(AR*lambda));
C_l_Beta_v        = -a_v*(S_v/S_w)*(z_v/b_w)*(eta_v);
C_l_Beta          = (C_l_Beta_w+C_l_Beta_w_Gamma0+C_l_Beta_v+delta_C_l_Beta_1+delta_C_l_Beta_2)*0.5; %Correction Factor
C_n_Beta_t        = -(l_v/b_w)*C_Y_Beta_t;
C_n_Beta          = C_n_Beta_t-K_n*(S_s/S_w)*(l_b/b_w);
C_l_p_w           = C_l_p_a_0w*((AR+4)/(4+AR*2*pi/a_0_w))-C_D/8;
C_l_p_t           = 0.5*(S_t/S_w)*((b_t/b_w)^2)*(C_l_p_a_0t)*((AR_t+4*cosd(Lambda_t))/(4*cosd(Lambda_t)+AR_t*2*pi/a_0_t));
C_l_p_v           = 2*(z_v/b_w)*(z_v/b_w)*C_Y_Beta_t;
C_l_p             = C_l_p_t+C_l_p_v+C_l_p_w;
C_Y_p             = (C_Y_p_over_C_L)*C_L + C_l_p*(Delta_C_Y_p_Gamma_over_C_l_p_Gamma0);
C_n_p_w           = C_L*((AR+4)/(AR+4*cosd(Lambda)))*(1+0.5*(1+cosd(Lambda)/AR)*(tand(Lambda))^2)*(C_n_p_over_C_L_Lambda0);
C_n_p_v           = a_v*(S_v/S_w)*(1/b_w)*(z_v*sind(alpha_v)+l_v*cosd(alpha_v))*((2/b_w)*(z_v*cosd(alpha_v)-l_v*sind(alpha_v))-(dsigma1+dsigma2));
C_n_p             = (C_n_p_w+C_n_p_v)*(0.75); %Correction Factor Added
C_Y_r_t           = -2*(l_v/b_w)*C_Y_Beta_t;
C_Y_r_w           = 0.143*C_L-0.05;
C_Y_r             = C_Y_r_t+C_Y_r_w;
C_l_r             = C_L/4-2*(l_v/b_w)*(z_v/b_w)*C_Y_Beta_t;
C_n_r_t           = -2*(l_v/b_w)*C_n_Beta_t;
C_n_r_w           = -(0.33*C_D_0*(1+3*lambda)/(2+2*lambda)+0.02*C_L*C_L*(1-(AR-6)/13-(1-lambda)/2.5));
C_n_r             = (C_n_r_t+C_n_r_w)*(1.5); %Correction Factor Added
C_Y_delta_A       = 0;
C_l_delta_A       = (C_L_delta_A_over_tauo-C_L_delta_A_over_taui)*tau1*(0.55); %Correction Factor Added
C_n_delta_A       = 2*K*C_L*C_l_delta_A*(0.5); %Correction Factor Added
C_Y_delta_R       = a_v*tau2*(S_v/S_w)*(0.9); %Correction Factor Added
C_l_delta_R       = a_v*tau2*(S_v/S_w)*(z_v/b_w);
C_n_delta_R       = -a_v*tau2*(S_v/S_w)*(l_v/b_w)*(eta_v)*(2); %Correction Factor Added
disp('Calculation Complete')  


% Problem #8
Calculated_Values = [C_L C_D C_m C_L_u C_D_u C_m_u C_L_alpha C_D_alpha C_m_alpha C_L_alphadot C_D_alphadot C_m_alphadot C_L_q C_D_q C_m_q C_L_delta_e C_D_delta_e C_m_delta_e C_Y_Beta C_l_Beta C_n_Beta C_Y_p C_l_p C_n_p C_Y_r C_l_r C_n_r C_Y_delta_A C_l_delta_A C_n_delta_A C_Y_delta_R C_l_delta_R C_n_delta_R].';
Stability_Derivative = ["CL" "CD" "CM" "CL_u" "CD_u" "CM_u" "CL_alpha" "CD_alpha" "CM_alpha" "CL_alpha_dot" "CD_alpha_dot" "CM_alpha_dot" "CL_q" "CD_q" "CM_q" "CL_deltaE" "CD_deltaE" "CM_deltaE" "CY_beta" "Cl_beta" "Cn_beta" "CY_p" "Cl_p" "Cn_p" "CY_r" "Cl_r" "Cn_r" "CY_deltaA" "Cl_deltaA" "Cn_deltaA" "CY_deltaR" "Cl_deltaR" "Cn_deltaR"].';
Published_Values = [0.445 0.0224 NaN NaN NaN NaN 4.8762 0.212 -1.5013 NaN NaN -4.10 NaN NaN -12.05 0.328 0 -0.971 -0.6532 -0.13752 0.12319 NaN -0.416 -0.0307 NaN 0.132 -0.161 0 0.08308 -0.00354 0.18651 0.019195 -0.08337].';
Percent_Difference = 100*abs((Calculated_Values-Published_Values)./Published_Values);
Tab = table(Stability_Derivative,Published_Values,Calculated_Values,Percent_Difference)


% Problem #5
%% Empirical Factor for fuselage or nacelle P.49 Figure 9 
function [k_f] = figure9(x_ac,l_b)
x0  = (x_ac/l_b)*100;
x   = [10 20 30 40 50 60]';
y   = [0.08 0.3 0.6 1.0 1.6 2.8]';
f   = fit(x,y,'poly3');
f   = coeffvalues(f);
k_f = f(1)*x0^3+f(2)*x0^2+f(3)*x0+f(4);
end 
%% Elevator Effectivness P.64 Figure 13 
function [dalpha_t_over_ddelta_e] = figure13(S_e,S_t)
x0 = S_e/S_t;
x  = [0 .1 .2 .3 .4 .5 .6 .7]';
y  = [0 0.27 0.41 0.51 0.6 0.67 0.73 0.78]';
f  = fit(x,y,'poly3');
f  = coeffvalues(f);
dalpha_t_over_ddelta_e = f(1)*x0^3+f(2)*x0^2+f(3)*x0+f(4);
end 
%% Flap Effectivness parameter P.65 Figure 14 
function [C_l_delta_e] = figure14(c_a,c,AR_t)
x0   = c_a/c;
x    = [0 .2 .4 .6 .8 1.0 ]';
y1   = [0 .031 .054 .077 .097 0.11]';
y2   = [0 .031 .053 .066 .071 .072]';
if AR_t > 0.5 % High AR 
   y = y1; 
else % Low AR 
   y = y2;
end 
f    = fit(x,y,'poly3');
f    = coeffvalues(f);
C_l_delta_e = f(1)*x0^3+f(2)*x0^2+f(3)*x0+f(4);
end 
%% Wing-Fuselage Interference Fator P.73 Figure 23
function [K_i] = figure23 (z_w,d)
x0 = z_w/(d/2);
x = [0 .2 .4 .6 .8]';
    if z_w<0
        x = -x;
        y = [1.0 1.175 1.35 1.52 1.68]';
    else
        y = [1.0 1.1 1.2 1.3 1.4]'; 
    end 
        f = fit(x,y,'poly3');
        f = coeffvalues(f);
        K_i = f(1)*x0^3+f(2)*x0^2+f(3)*x0+f(4);
end 
%% K = f(ratio of vertical tail span to fuselage diameter around tail) P.74 Figure 24 
function [k] = figure24 (b_v,r_1)
x0 = b_v/(2*r_1);
    if x0 <= 2
        k = 0.77;
    elseif x0 >= 3.5
        k = 1;
    else 
        f = fit([2;3.5],[.77;1],'poly1');
        f = coeffvalues(f);
        k = f(1)*x0+f(2);
    end 
end
%% Cl_Beta_over_Gamma for various AR and lambda P.77 Figure 26
function [C_l_Beta_over_Gamma] = figure26 (AR,lambda)
x0      = AR;
lambdas = [1 0.5 0];
n       = abs(lambda-lambdas)==min(abs(lambda-lambdas));
x       = [0 2 4 6 8]';
y(:,1)  = [0 -.000095 -.000175 -.000215 -.000255]';
y(:,2)  = [0 -.000095 -.000170 -.000208 -.000240]';
y(:,3)  = [0 -.000090 -.000140 -.000170 -.000192]';
y       = y(:,n);
f       = fit(x,y,'poly3');
f       = coeffvalues(f);
C_l_Beta_over_Gamma = f(1)*x0^3+f(2)*x0^2+f(3)*x0+f(4);
end 
%% a_v as function of vertical tail AR P.81 Figure 28
function [a_v] = figure28(Ae)
x0      = Ae;
x       = [0 1 2 3 4 5 6]';
y       = [0 0.03 0.045 0.057 0.064 0.07 0.073]';
f       = fit(x,y,'poly2');
f       = coeffvalues(f);
a_v     = f(1)*x0^2+f(2)*x0^1+f(3);
end 
%% Delta_C_Y_p_Gamma_over_C_l_p_Gamma0 as function of dihedral P.87 Figure 33
function [Delta_C_Y_p_Gamma_over_C_l_p_Gamma0] = figure33(Gamma)
x0 = Gamma;
x = [-10 -5 0 5 10]';
y = [-0.65 -0.38 -0.08 0.22 0.52]';
f = fit(x,y,'poly1');
f = coeffvalues(f);
Delta_C_Y_p_Gamma_over_C_l_p_Gamma0  = f(1)*x0+f(2);
end 
%% Wing Contribution to C_l_p for 0 or small sweep Lambda P.88 Figure 34 
function [C_l_p_over_C_L0] = figure34(Lambda,AR)
Lambda  = Lambda*180/pi;
x0      = AR; 
Lambdas = [0 0.5 1];
n       = abs(Lambda-Lambdas)==min(abs(Lambda-Lambdas));
x       = [0 2 4 6 8 10]';
y(:,1)  = [0 -.16 -.25 -.31 -.35 -.36]';
y(:,2)  = [0 -.19 -.315 -.405 -.48 -.52]';
y(:,3)  = [0 -.19 -.32 -.44 -.515 -.59]';
y       = y(:,n);
f       = fit(x,y,'poly3');
f       = coeffvalues(f);
C_l_p_over_C_L0 = f(1)*x0^3+f(2)*x0^2+f(3)*x0+f(4);
end 
%% Wing Contribution to C_n_p/CL for zero sweep P. 91 Figure 35 
function [C_n_p_over_C_L_Lambda0] = figure35(Lambda,AR)
Lambda  = Lambda*180/pi;
x0      = AR; 
Lambdas = [1 0.5];
n       = abs(Lambda-Lambdas)==min(abs(Lambda-Lambdas));
x       = [2 4 6 8 10 12 14 16]';
y(:,1)  = [-.02 -.039 -.05 -.06 -.067 -.075 -.082 -.087]';
y(:,2)  = [-.02 -.0395 -.052 -.063 -.07 -.078 -.082 -.085]';
y       = y(:,n);
f       = fit(x,y,'poly3');
f       = coeffvalues(f);
C_n_p_over_C_L_Lambda0 = f(1)*x0^3+f(2)*x0^2+f(3)*x0+f(4);
end 
%% Effect of wing sidewash on vertical tail P. 92 Figure 36
function [dsigma_1_over_dpb2U] = figure36(Lambda,z_v,b_w)
% Lambda in degree input
Lambda  = Lambda*180/pi;
x0      = z_v/(b_w/2); 
Lambdas = [1 0.5 0];
n       = abs(Lambda-Lambdas)==min(abs(Lambda-Lambdas));
x       = [.6 .55 .5 .45 .4 .35 .3 .25 .2 .15 .1 .05 0]';
y(:,1)  = [.1 .11 .12 .135 .145 .160 .175 .2 .22 .248 .275 .27 0]';
y(:,2)  = [.09 .11 .11 .135 .145 .160 .175 .2 .22 .248 .275 .31 0]';
y(:,3)  = [.075 .085 0.1 .115 .125 .150 .17 .195 .22 .245 .295 .345 0]';
y       = y(:,n);
f       = fit(x,y,'poly3');
f       = coeffvalues(f);
dsigma_1_over_dpb2U = f(1)*x0^3+f(2)*x0^2+f(3)*x0+f(4);
end 
%% tau as function of c_a to c_w ratio P.102 Figure 41
function [tau1] = figure41(c_a,c)
x0   = c_a/c; 
x    = [0.1 1.5 0.2 0.25 0.3]';
y    = [0.235 0.315 0.375 0.435 0.475]';
f    = fit(x,y,'poly3');
f    = coeffvalues(f);
tau1 = f(1)*x0^3+f(2)*x0^2+f(3)*x0+f(4);
end
%% tau as function of c_a to c_w ratio P.102 Figure 42
function [C_L_delta_A_over_tau] = figure42(lambda,Y_i,b_a,b_w)
Lambda  = lambda*180/pi;
Lambdas = [1 0.5 0];
n       = abs(Lambda-Lambdas)==min(abs(Lambda-Lambdas));
x       = [0 0.1 0.2 0.3 0.4 0.6 0.8 0.9 1.0]';
% AR = 6 
x0      = (Y_i+b_a)/(b_w/2);
y(:,1)  = [0 0.02 0.05 0.09 0.15 0.37 0.6 0.7 0.79]';
y(:,2)  = [0 0.02 0.05 0.09 0.15 0.37 0.59 0.68 0.76]';
y(:,3)  = [0 0.02 0.05 0.09 0.15 0.37 0.58 0.65 0.72]';
y       = y(:,n);
f       = fit(x,y,'poly4');
f       = coeffvalues(f);
C_L_delta_A_over_tau = f(1)*x0^4+f(2)*x0^3+f(3)*x0^2+f(4)*x0+f(5);
end 
%% Empirical Factor K of eta for taper ratio = 0.5 Figure 43
function [K] = figure43(Y_i,b_w,AR)
x0     = Y_i/(b_w/2);
ARs    = [4 6 8];
n      = abs(AR-ARs)==min(abs(AR-ARs));
x      = [0 0.25 0.5 0.75]';
y(:,1) = [-0.24 -0.245 -0.25 -0.26]';
y(:,2) = [-0.18 -0.18 -0.185 -0.195]';
y(:,3) = [-0.13 -0.13 -0.135 -0.15]';
y      = y(:,n);
f      = fit(x,y,'poly3');
f      = coeffvalues(f);
K      = f(1)*x0^3+f(2)*x0^2+f(3)*x0+f(4);
end 
%% tau as a function of rudder area to vertical tail area P.107 Figure 45/46
function [tau2] = figure45(S_r,S_v)
x0   = S_r/S_v;
x    = [0 .1 .2 .3 .4 .5 .6 .7]';
y    = [0 .3 .4 .5 .6 .68 .63 .78]';
f    = fit(x,y,'poly3');
f    = coeffvalues(f);
tau2 = f(1)*x0^3+f(2)*x0^2+f(3)*x0+f(4);
end 
% \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\