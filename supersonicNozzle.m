% Supersonic Nozzle (Method of Characteristics)
clear all, close all, clc

%%PROMPT
try
    prompt1 = 'Mach Number?: ';
    Mexit = input(prompt1);
catch
    warning('Invalid Mach Number. Assigning a value of 2.5.');
    Mexit = 2.5;
end

try
    prompt2 = 'Number of characteristic lines: ';
    n = input(prompt2);
catch
    warning('Invalid characteristic lines. Assigning a value of 7');
    n = 7;
end

prompt3 = 'Plot the characteristic lines? [Y/N]: ';
plotLines = input(prompt3,'s');
lines = 0;
if (strncmpi(plotLines,'y',1) | strncmpi(plotLines,'Y',1))
    lines = 1;
end

prompt4 = 'Plot the gradients? [Y/N]: ';
plotGradients = input(prompt4,'s');
gradients = 0;
if (strncmpi(plotGradients,'y',1) | strncmpi(plotGradients,'Y',1))
    gradients = 1;
end

prompt5 = 'Output full tables, partial tables, or no tables? [F/P/N]: ';
outputTables = input(prompt5,'s');
tables = 0;
if (strncmpi(outputTables,'f',1) | strncmpi(outputTables,'F',1))
    tables = 2;
elseif (strncmpi(outputTables,'p',1) | strncmpi(outputTables,'P',1))
    tables = 1;
end

%% Inputs and Function Definitions

l = 1; % Arbitrary half throat width
% Define necessary inline functions; "_" at the end of an expression
    % indicates a function rather than a known variable

Kneg_ = @(Theta,nu) Theta + nu;
Kpos_ = @(Theta,nu) Theta - nu;
Theta_ = @(Kneg,Kpos) .5*(Kneg+Kpos);
Nu_ = @(Kneg,Kpos) .5*(Kneg-Kpos);
T_rat_ = @(M) (1+.2.*M.^2).^-1;
p_rat_ = @(M) (1+.2.*M.^2).^(-1.4/.4);
rho_rat_ = @(M) (1+.2.*M.^2).^(-1/.4);
% Output/Input nu values are in degrees
nu_ = @(M) (sqrt(6).*atand(sqrt((M.^2-1)/6))-atand(sqrt(M.^2-1)));
A = 1.3604;
B = .0962;
C = -.5127;
D = -.6722;
E = -.3278;
numax = pi/2*(sqrt(6)-1)*180/pi;
M_ = @(nu) (1+A*(nu./numax).^(2/3)+B*(nu./numax).^(4/3)+C*(nu./numax).^(2))./...
    (1+D*(nu./numax).^(2/3)+E*(nu./numax).^(4/3));

%% Table of Values for All points

% 1) Table of values with has 11 columns; 
    % Columns 1 - 4 correspond to 
        % Kneg, Kpos, Theta, and nu in that order, where 2 are needed at
        % previous points to calculate those at the next;
    % Columns 5 and 6 correspond to
        % M and mu calculated from nu
    % Columns 7 - 9 correspond to
        % Ratios of T, p, and rho to their stagnation values
    % Columns 10 and 11 are the x and y coor. of the points;
        % to be calculated later

% 2) Table has R = (n+1)(n+2)/2-1 rows because there are R point
    % where the characteristic lines intersect a wall, the center line,
    % or other characteristic lines

Points = (n+1)*(n+2)/2-1;
Data = zeros(Points,11);

% 3) The maximum theta is equal to half of the nu_(Mexit)
% 4) Theta(1) to Theta(n) are equally spaced between 0 and the Theta_max
    % e.g. Theta(n) = Theta_max, and Theta(1) = Theta_max/n

Theta_max = nu_(Mexit)/2;
Data(1:n,3) = linspace(Theta_max/n,Theta_max,n);
Data(1:n,4) = Data(1:n,3);

Pcount1 = 1;

for j = 1:n  
    Data(j,1) = Kneg_(Data(j,3),Data(j,4));
    Data(j,2) = Kpos_(Data(j,3),Data(j,4));
    Pcount1 = Pcount1 + 1; 
end

Data(n+1,1:4) = Data(n,1:4);

for k = 1:(n-1);
    Pcount1 = Pcount1 + 1;
    Data(Pcount1,3) = 0;
    Data(Pcount1,1) = Data(Pcount1-n+k-1,1);
    Data(Pcount1,2) = -Data(Pcount1,1);
    Data(Pcount1,4) = Data(Pcount1,1);
    for j = (Pcount1+1):(Pcount1+n-k-1)
        Pcount1 = Pcount1 + 1;
        Data(Pcount1,1) = Data(Pcount1-n+k-1,1);
        Data(Pcount1,2) = Data(Pcount1-1,2);
        Data(Pcount1,3) = Theta_(Data(Pcount1,1),Data(Pcount1,2));
        Data(Pcount1,4) = Nu_(Data(Pcount1,1),Data(Pcount1,2));
    end
    Pcount1 = Pcount1 + 1;
    Data(Pcount1,1:4) = Data(Pcount1-1,1:4);
end

Data(:,5) = M_(Data(:,4));
Data(:,6) = asind(1./Data(:,5));
Data(:,7) = T_rat_(Data(:,5));
Data(:,8) = p_rat_(Data(:,5));
Data(:,9) = rho_rat_(Data(:,5));
Data;

%% Find the Coordinates for Each Point

% Define formulas for slope
slpneg = @(a,b) tand(.5*(Data(a,3)+Data(b,3)) - .5*(Data(a,6)+Data(b,6)));
slppos = @(a,b) tand(.5*(Data(a,3)+Data(b,3)) + .5*(Data(a,6)+Data(b,6)));
slpwall = @(a,b) tand(.5*(Data(a,3)+Data(b,3)));

Pcount2 = 1;

% For the Points along the first Kpos line
Data(1,10) = -l/tand(Data(1,3)-Data(1,6));
for j = 2:n
    Pcount2 = Pcount2 + 1;
    Data(j,10) = (l-Data(j-1,11)+slppos(j,j-1)*Data(j-1,10))...
        /(slppos(j,j-1)-tand(Data(j,3)-Data(j,6)));
    Data(j,11) = Data(j-1,11)+slppos(j,j-1)*(Data(j,10)...
        -Data(j-1,10));
end
Pcount2 = Pcount2 + 1;

Data(Pcount2,10) = (l-Data(Pcount2-1,11)+slppos(Pcount2,Pcount2-1)*...
    Data(Pcount2-1,10))/(slppos(Pcount2,Pcount2-1)-tand(Theta_max));
Data(Pcount2,11) = l + tand(Theta_max)*Data(Pcount2,10);

for k = 1:(n-1);
    Pcount2 = Pcount2 + 1;
    Data(Pcount2,10) = Data(Pcount2-n+k-1,10)-Data(Pcount2-n+k-1,11)/...
        slpneg(Pcount2,Pcount2-n+k-1);
    
    for j = (Pcount2+1):(Pcount2+n-k-1)         
        Pcount2 = Pcount2 + 1;
%         Data(Pcount2,1) = Data(Pcount2-n+k-1,1);
%         Data(Pcount2,2) = Data(Pcount2-1,2);  
        Data(j,10) = (Data(j-n+k-1,11)-Data(j-1,11)+...
            slppos(j,j-1)*Data(j-1,10)-...
            slpneg(j,j-n+k-1)*Data(j-n+k-1,10))/...
            (slppos(j,j-1)-slpneg(j,j-n+k-1));
        Data(j,11) = Data(j-1,11)+slppos(j,j-1)...
            *(Data(j,10)-Data(j-1,10));  
     end
     Pcount2 = Pcount2 + 1;
     Data(Pcount2,10) = (Data(Pcount2-n+k-1,11)-Data(Pcount2-1,11)+...
            slppos(Pcount2,Pcount2-1)*Data(Pcount2-1,10)-...
            slpwall(Pcount2,Pcount2-n+k-1)*Data(Pcount2-n+k-1,10))/...
            (slppos(Pcount2,Pcount2-1)-slpwall(Pcount2,Pcount2-n+k-1));
        Data(Pcount2,11) = Data(Pcount2-1,11)+slppos(Pcount2,Pcount2-1)...
            *(Data(Pcount2,10)-Data(Pcount2-1,10));
end

Rowlabel = (1:Pcount1);
Light_output = [Data(:,5),Data(:,7),Data(:,8),Data(:,9),Data(:,10),Data(:,11)];
if tables == 1
    printmat(Light_output, 'Partial Output', num2str(Rowlabel) , ...
    'M T/T_o p/p_o rho/rho_o x-coor y-coor' )
end
if tables == 2
    printmat(Data, 'Full Output', num2str(Rowlabel) , ...
    'K- K+ Theta nu M mu T/T_o p/p_o rho/rho_o x-coor y-coor' )
end

%% Flow Field Mapping
% Define extra points for whole flow field with assumed M at those points
    % Need 1 point at throat center, 1 point at throat edge, 1 point
    % halfway between throat edge and pt. 1, 1 pt. at the centerline exit,
    % and a few points along the final Kpos line, all with appropriate data

    if gradients == 1
        DataXtra = zeros(10,11);
        DataXtra(2,11) = 1;
        DataXtra(3,10:11) = [Data(1,10)/2,l/2];
        DataXtra(1:3,5) = 1;
        DataXtra(4,10:11) = [Data(end,10),0];
        DataXtra(5:10,10:11) = [linspace(Data(end-1,10),Data(end,10),6)',...
            linspace(Data(end-1,11),Data(end,11),6)'];
        DataXtra(4:10,5) = Data(end,5);
        
        %Define other properties based on M at each point
        DataXtra(:,7) = T_rat_(DataXtra(:,5));
        DataXtra(:,8) = p_rat_(DataXtra(:,5));
        DataXtra(:,9) = rho_rat_(DataXtra(:,5));
        
        warning('off','all')
        
        xint = [Data(:,10);DataXtra(:,10)];
        yint = [Data(:,11);DataXtra(:,11)];
        Mint = [Data(:,5);DataXtra(:,5)];
        Tint = [Data(:,7);DataXtra(:,7)];
        pint = [Data(:,8);DataXtra(:,8)];
        rhoint = [Data(:,9);DataXtra(:,9)];
        
        deltax = linspace(0,Data(Pcount1,10),500);
        deltay = linspace(0,Data(Pcount1,11),500);
        
        Mscint = scatteredInterpolant(xint,yint,Mint,'linear','none');
        Tscint = scatteredInterpolant(xint,yint,Tint,'linear','none');
        pscint = scatteredInterpolant(xint,yint,pint,'linear','none');
        rhoscint = scatteredInterpolant(xint,yint,rhoint,'linear','none');
        
        [xgrid,ygrid] = meshgrid(deltax,deltay);
        
        figure('units','normalized','outerposition',[0 0 1 .5])
        colormap jet
        subplot(2,2,1)
        hold on
        mesh(xgrid,ygrid,Mscint(xgrid,ygrid))
        mesh(xgrid,-ygrid,Mscint(xgrid,ygrid))
        
        view(2)
        xlabel('x position')
        ylabel('y position')
        title(['Mach Flow Field for ' num2str(n) ' Char. Lines and M_e = '...
            num2str(Mexit)])
        set(gca,'fontsize',14)
        colorbar
        
        subplot(2,2,2)
        hold on
        mesh(xgrid,ygrid,Tscint(xgrid,ygrid))
        mesh(xgrid,-ygrid,Tscint(xgrid,ygrid))
        
        view(2)
        xlabel('x position')
        ylabel('y position')
        title(['T/T_0 Flow Field for ' num2str(n) ' Char. Lines and M_e = '...
            num2str(Mexit)])
        set(gca,'fontsize',14)
        colorbar
        
        subplot(2,2,3)
        hold on
        mesh(xgrid,ygrid,pscint(xgrid,ygrid))
        mesh(xgrid,-ygrid,pscint(xgrid,ygrid))
        
        view(2)
        view(2)
        xlabel('x position')
        ylabel('y position')
        title(['p/p_0 Flow Field for ' num2str(n) ' Char. Lines and M_e = '...
            num2str(Mexit)])
        set(gca,'fontsize',14)
        colorbar
        
        subplot(2,2,4)
        hold on
        mesh(xgrid,ygrid,rhoscint(xgrid,ygrid))
        mesh(xgrid,-ygrid,rhoscint(xgrid,ygrid))
        
        view(2)
        xlabel('x position')
        ylabel('y position')
        title(['\rho/\rho_0 Flow Field for ' num2str(n) ' Char. Lines and M_e = '...
            num2str(Mexit)])
        set(gca,'fontsize',14)
        colorbar
    end

%% Distinguishing points along the wall, centerline points, 
    % and free points along K- and K+ lines

% Wall points
wallpoints = zeros(n,1);
increment = n;
wallpoints(1,1) = n+1;
for i = 2:n;
    wallpoints(i,1) = wallpoints(i-1,1) + increment;
    increment = increment - 1;
end

% Center points
centerpoints = zeros(n,1);
centerpoints(1,1) = 1;
for i = 2:n;
    centerpoints(i,1) = wallpoints(i-1,1) + 1;
end

% K+ points
kpospoints = zeros(n+1);
for i = 1:n
    increment = 0;
    for j = 1:n+1
        if (increment + centerpoints(i,1) <= wallpoints(i,1))
            kpospoints(i,j) = centerpoints(i,1) + increment;
        end
        increment = increment + 1;
    end
end

% K- points
knegpoints = zeros(n) ;
maxCount = 2;
for i = 1:n
    knegpoints(i,1) = i;
    for j = 2:maxCount
        %knegpoints(i,j) = centerpoints(i,1)
        knegpoints(i,j) = knegpoints(i,j-1) + (n - j + 2);
    end
    maxCount = maxCount + 1;
end
for i = 1:n
    knegpoints(i,i+1) = 0;
end

%% Plotting the Lines
if lines == 1
    figure('units','normalized','outerposition',[0 .53 1 .47])
    hold on
    
    for i = 1:n
        line = kpospoints(i,:);
        line(line==0)=[];
        x = Data(line,10);
        y = Data(line,11);
        plot(x,y,'-g');
        plot(x,-y,'-g');
        
        line = knegpoints(i,:);
        line(line==0)=[];
        x = [0; Data(line,10)];
        y = [1; Data(line,11)];
        plot(x,y,'-g');
        plot(x,-y,'-g');
    end
    
    x = [0; Data(wallpoints(:,1),10)];
    y = [1; Data(wallpoints(:,1),11)];
    plot(x,y,'-k');
    plot(x,-y,'-k')
    
    plot([0,Data(end,10)],[0,0],'-.k')
    
    xlabel('x position')
    ylabel('y position')
    title(['' num2str(n) ' Char. Lines for M_e = '...
        num2str(Mexit)])
    set(gca,'fontsize',14)
    axis equal
    xlim([0 7])
    ylim([-2.5 2.5])
    grid on; grid minor;
end
