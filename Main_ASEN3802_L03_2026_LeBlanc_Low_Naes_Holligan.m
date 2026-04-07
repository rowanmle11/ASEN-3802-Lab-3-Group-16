% Contributors: Rowan LeBlanc

clear;clc;close all

%% Task 1 NACA 4-Digit Airfoil Generator

c=1;
N=100;

% NACA 0021
[x1,y1,xc_0021,yc_0021] = NACA_Airfoils(0,0,0.21,c,N);

% NACA 2421
[x2,y2,xcam,ycam] = NACA_Airfoils(0.02,0.4,0.21,c,N);


figure
hold on
grid on
plot(x1, y1, 'b-')
plot(x2, y2, 'r-')
plot(xcam,ycam,'k-')
axis equal
xlabel('x/c')
ylabel('y/c')
legend('NACA 0021','NACA 2421','Camber Line')
title('NACA Airfoil Shape')

%% Task 2 Convergence Study on Number of Panels Required

alpha = 12;

Num_Panels = 500;
chord = 1;
cl = zeros(1,Num_Panels);
    count = 0;

    for i = 2:Num_Panels+1
        count = count + 1
        [x,y,xc,yc] = NACA_Airfoils(0,0,0.12,chord,i);
        cl(i-1) = Vortex_Panel(x,y,alpha);
    end

c_intersect = (cl(length(cl))- 0.01*cl(length(cl))) * ones(1,Num_Panels);
x_axis = 2:Num_Panels+1;
[xid,yid] = polyxpoly(x_axis,cl,x_axis,c_intersect);

disp('Number of Panels to 1% error:')
disp(ceil(xid))

figure();
hold on
plot(x_axis,cl, 'LineWidth',2);
xlabel('Number of Panels')
ylabel('Cl Values')
yline(cl(length(cl)) + 0.01*cl(length(cl)), 'r--')
yline(cl(length(cl)) - 0.01*cl(length(cl)), 'r--')
xline(ceil(xid),'k', 'LineWidth', 1.5);
plot(ceil(xid),yid,'Marker','o','MarkerSize',8,'LineWidth',2,'Color','b');
xlim([0,Num_Panels]);
legend('C_{L}', '+1% C_{L final}','-1% C{L final}','Min Number of Panels = 29', '1% intersect','Interpreter', 'tex', 'Location','Best')
title('Cl vs. Number of Panels')
hold off

cl_exact = cl(end);
panels_exact = x_axis(end);

min_panels_pred = ceil(xid(1));
cl_pred = interp1(x_axis, cl, min_panels_pred);

rel_error = abs(cl_pred - cl_exact)/abs(cl_exact);

Results = table(cl_exact, panels_exact, cl_pred, rel_error, min_panels_pred,'VariableNames', {'C_{L Exact}','Panels_Exact','C{L Predicted}','Relative_Error','Min_Panels_Predicted'});

disp(Results)

%% Task 3 Effect of Airfoil Thickness on Lift

alpha = -5:1:10;
N = 50;
c = 1;

airfoils = {'NACA0006','NACA0012','NACA0018'};
t = [0.06,0.12,0.18];
m = 0;
p = 0;

cl_panel = zeros(length(t),length(alpha));

for i=1:length(t)
    [xb,yb,~,~] = NACA_Airfoils(m,p,t(i),c,N);
    for j=1:length(alpha)
        cl_panel(i,j) = Vortex_Panel(xb,yb,alpha(j));
    end
end

% Thin Airfoil Theory
a0 = 2*pi*(pi/180);
cl_tat = a0*alpha;

% Experimental Data

% Plots

figure
hold on
grid on
plot(alpha,cl_tat,'k--','DisplayName','Thin Airfoil Theory')

plot(alpha,cl_panel(1,:),'b-','DisplayName','NACA 0006')
plot(alpha,cl_panel(2,:),'r-','DisplayName','NACA 0012')
plot(alpha,cl_panel(3,:),'g-','DisplayName','NACA 0018')

% experimental data plots

xlabel('Angle of Attack, \alpha (degrees)')
ylabel('Sectional Lift Coefficient, c_l')
title('Effect of Airfoil Thcikness on Lift Coefficient')
legend('Location','best')
hold off

linear_idx = find(alpha >= -5 & alpha <= 5);

fprintf('\nZero-Lift AoA (degrees):\n')
fprintf('%-15s | %-15s | %-15s | %-15s\n','Airfoil','TAT','Vortex Panel','Experimental')
fprintf(repmat('-',1,65),'\n')

fprintf('\nLift Slope (1/degrees):\n')
fprintf('%-15s | %-15s | %-15s | %-15s\n','Airfoil','TAT','Vortex Panel','Experimental')
fprintf(repmat('-',1,65),'\n')

for i = 1:length(t)
    P_panel = polyfit(alpha(linear_idx),cl_panel(i,linear_idx),1);
    slope_panel = P_panel(1);
    alpha_L0_panel = -P_panel(2)/P_panel(1);

    fprintf('\n%-15s | %-15.2f | %-15.2f | %-15s\n',airfoils{i},0,alpha_L0_panel,' ')
end

%% Task 4 Effect of Airfoil Camber on Lift

N = ceil(xid);

alpha_vals = -100:100;
chord = 1;
cl1 = zeros(1,length(alpha_vals));
cl2 = zeros(1,length(alpha_vals));
cl3 = zeros(1,length(alpha_vals));

for i = 1:length(alpha_vals)

    [x1,y1,xc1,yc1] = NACA_Airfoils(0,0,0.12,chord,N);
    cl1(i) = Vortex_Panel(x1,y1,alpha_vals(i));

    [x2,y2,xc2,yc2] = NACA_Airfoils(0.02,0.4,0.12,chord,N);
    cl2(i) = Vortex_Panel(x2,y2,alpha_vals(i));
    
    [x3,y3,xc3,yc3] = NACA_Airfoils(0.04,0.4,0.12,chord,N);
    cl3(i) = Vortex_Panel(x3,y3,alpha_vals(i));

end

c_intersect = zeros(1,length(alpha_vals));
x_axis = 1:length(alpha_vals);
[xid1,yid1] = polyxpoly(x_axis,cl1,x_axis,c_intersect);
[xid2,yid2] = polyxpoly(x_axis,cl2,x_axis,c_intersect);
[xid3,yid3] = polyxpoly(x_axis,cl3,x_axis,c_intersect);

zero_lift_AOA1 = cl1(ceil(xid1))
zero_lift_AOA2 = cl2(ceil(xid2))
zero_lift_AOA3 = cl3(ceil(xid3))

figure();
hold on
plot(alpha_vals,cl1, 'Color', 'b', 'LineWidth', 2)
plot(alpha_vals,cl2, 'Color', 'r', 'LineWidth', 2)
plot(alpha_vals,cl3, 'Color', 'magenta', 'LineWidth', 2)
xlabel('AOA \alpha [Deg]', 'Interpreter', 'tex')
ylabel('C_{L}')
legend('NACA 0012', 'NACA 2412', 'NACA 4412', 'Location', 'Best')
title('C_{L} vs. \alpha', 'Interpreter', 'tex')

hold off

cl1_id = find(cl1 == 0, 1);
cl2_id = find(cl2 == 0, 1);
cl3_id = find(cl3 == 0, 1);

%% FUNCTIONS

function [x_b, y_b, x, y_c] = NACA_Airfoils(m,p,t,c,N)

    % cosine (equiangular spacing)
    beta = linspace(0,pi,N+1);
    x = (c/2)*(1-cos(beta));

    y_t = (t/0.2)*c*(0.2969*sqrt(x/c)-0.1260*(x/c)-0.3516*(x/c).^2+0.2843*(x/c).^3-0.1036*(x/c).^4);

    y_c = zeros(size(x));
    dyc_dx = zeros(size(x));

    forward = x < p*c;

    if p > 0 && p < 1
        y_c(forward) = (m/p^2)*(2*p*(x(forward)/c) - (x(forward)/c).^2) * c;
        dyc_dx(forward) = (2*m / p^2) * (p-x(forward)/c);

        aft = ~forward;
        y_c(aft) = (m/(1-p)^2) * ((1-2*p) + 2*p*(x(aft)/c)-(x(aft)/c).^2)*c;
        dyc_dx(aft) = (2*m / (1-p)^2)*(p-x(aft)/c);
    end

    xi = atan(dyc_dx);

    % upper
    xU = x-y_t.*sin(xi);
    yU = y_c+y_t.*cos(xi);
    % lower
    xL = x+y_t.*sin(xi);
    yL = y_c-y_t.*cos(xi);

    x_lower = xL(end:-1:1);
    y_lower = yL(end:-1:1);

    x_upper = xU(2:end);
    y_upper = yU(2:end);

    x_b = [x_lower,x_upper];
    y_b = [y_lower,y_upper];
end

function [CL] = Vortex_Panel(XB,YB,ALPHA)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input: %
% %
% XB = Boundary Points x-location %
% YB = Boundary Points y-location %
% VINF = Free-stream Flow Speed %
% ALPHA = AOA %
% %
% Output: %
% %
% CL = Sectional Lift Coefficient %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%
% Convert to Radians %
%%%%%%%%%%%%%%%%%%%%%%

ALPHA = ALPHA*pi/180;

%%%%%%%%%%%%%%%%%%%%%
% Compute the Chord %
%%%%%%%%%%%%%%%%%%%%%

CHORD = max(XB)-min(XB);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine the Number of Panels %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M = max(size(XB,1),size(XB,2))-1;
MP1 = M+1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Intra-Panel Relationships: %
% %
% Determine the Control Points, Panel Sizes, and Panel Angles %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for I = 1:M
        IP1 = I+1;
        X(I) = 0.5*(XB(I)+XB(IP1));
        Y(I) = 0.5*(YB(I)+YB(IP1));
        S(I) = sqrt( (XB(IP1)-XB(I))^2 +( YB(IP1)-YB(I))^2 );
        THETA(I) = atan2( YB(IP1)-YB(I), XB(IP1)-XB(I) );
        SINE(I) = sin( THETA(I) );
        COSINE(I) = cos( THETA(I) );
        RHS(I) = sin( THETA(I)-ALPHA );
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inter-Panel Relationships: %
% %
% Determine the Integrals between Panels %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for I = 1:M
        for J = 1:M
            if I == J
            CN1(I,J) = -1.0;
            CN2(I,J) = 1.0;
            CT1(I,J) = 0.5*pi;
            CT2(I,J) = 0.5*pi;
            
            else
           
                A = -(X(I)-XB(J))*COSINE(J) - (Y(I)-YB(J))*SINE(J);
                B = (X(I)-XB(J))^2 + (Y(I)-YB(J))^2;
                C = sin( THETA(I)-THETA(J) );
                D = cos( THETA(I)-THETA(J) );
                E = (X(I)-XB(J))*SINE(J) - (Y(I)-YB(J))*COSINE(J);
                F = log( 1.0 + S(J)*(S(J)+2*A)/B );
                G = atan2( E*S(J), B+A*S(J) );
                P = (X(I)-XB(J)) * sin( THETA(I) - 2*THETA(J) ) ...
                + (Y(I)-YB(J)) * cos( THETA(I) - 2*THETA(J) );
                Q = (X(I)-XB(J)) * cos( THETA(I) - 2*THETA(J) ) ...
                - (Y(I)-YB(J)) * sin( THETA(I) - 2*THETA(J) );
                CN2(I,J) = D + 0.5*Q*F/S(J) - (A*C+D*E)*G/S(J);
                CN1(I,J) = 0.5*D*F + C*G - CN2(I,J);
                CT2(I,J) = C + 0.5*P*F/S(J) + (A*D-C*E)*G/S(J);
                CT1(I,J) = 0.5*C*F - D*G - CT2(I,J);
            end
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inter-Panel Relationships: %
% %
% Determine the Influence Coefficients %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for I = 1:M
        AN(I,1) = CN1(I,1);
        AN(I,MP1) = CN2(I,M);
        AT(I,1) = CT1(I,1);
        AT(I,MP1) = CT2(I,M);
        
            for J = 2:M
                AN(I,J) = CN1(I,J) + CN2(I,J-1);
                AT(I,J) = CT1(I,J) + CT2(I,J-1);
            end
    end

AN(MP1,1) = 1.0;
AN(MP1,MP1) = 1.0;

    for J = 2:M
        AN(MP1,J) = 0.0;
    end

RHS(MP1) = 0.0;

%%%%%%%%%%%%%%%%%%%%%%%%
% Solve for the gammas %
%%%%%%%%%%%%%%%%%%%%%%%%

GAMA = AN\RHS';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve for Tangential Veloity and Coefficient of Pressure %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for I = 1:M
        V(I) = cos( THETA(I)-ALPHA );
        
        for J = 1:MP1
            V(I) = V(I) + AT(I,J)*GAMA(J);
        end
        
        CP(I) = 1.0 - V(I)^2;
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve for Sectional Coefficient of Lift %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CIRCULATION = sum(S.*V);
CL = 2*CIRCULATION;

end
