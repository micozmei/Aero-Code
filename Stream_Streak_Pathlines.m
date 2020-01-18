clear all; close all; clc;

index = 0;
runtime = 3;
for t = 0:runtime/3:runtime
    index = index+1;
    [x,y] = meshgrid(0:0.1:1,0:0.1:1);
    % Velocity Profile
    u = x*t;
    v = -y;
    % Streamlines: dx/u = dy/v
    figure(1)
    subplot(2,2,index)
    % quiver(x,y,u,v)
    title(['Streamlines at t = ', num2str(t)])
    xlabel('x'), ylabel('y'), axis([0 1 0 1])
    startx = 0:0.05:1;
    starty = ones(size(startx));
    streamline(x,y,u,v,startx,starty)
end

% Extract Pathline Data for range of T_start
figure(2)
[P0x, P0y] = pathline(0,runtime);
P0 = plot(P0x, P0y,'k'); hold on
[P1x, P1y] = pathline(1,runtime);
P1 = plot(P1x, P1y,'b'); hold on
[P2x, P2y] = pathline(2,runtime);
P2 = plot(P2x, P2y,'g'); hold on
[P3x, P3y] = pathline(3,runtime);
P3 = plot(P3x, P3y,'r'); hold on

% Plot Streakline
xx = [P0x(61),P1x(41),P2x(21),P3x(1)];
xxD = double(xx);
yy = [P0y(61),P1y(41),P2y(21),P3y(1)];
yyD = double(yy);
z = 0:1:50;
p = interp1(xxD,yyD,z,'makima');
Strk = plot(z,p,'--m'); hold on
plot(P0x(61), P0y(61),'*k', P1x(41), P1y(41), '*b', P2x(21), P2y(21), '*g', P3x(1), P3y(1), '*r')
title('Pathlines and Streakline'), xlabel('x'), ylabel('y'), axis([0.5 46 0 0.5])
legend([P0, P1, P2, P3, Strk],{'t=0 Pathline','t=1 Pathline','t=2 Pathline','t=3 Pathline','Streakline'})

function [X,Y] = pathline(T_start,runtime)
    % Pathline: dx/dt = u, dy/dt = v
    timestep = T_start:0.05:T_start+runtime;
    syms x(t) y(t)
    ode1 = diff(x) == x*t;
    ode2 = diff(y) == -y;
    odes = [ode1; ode2];
    cond1 = x(T_start) == 0.5;
    cond2 = y(T_start) == 0.5;
    conds = [cond1; cond2];
    % Pathline Equation
    [xSol(t), ySol(t)] = dsolve(odes,conds);
    X = subs(xSol(t), t, timestep);
    Y = subs(ySol(t), t, timestep);
end
