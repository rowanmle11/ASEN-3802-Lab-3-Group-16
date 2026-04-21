% Contributors: Rowan LeBlanc, Landon Holligan, Erik Low, Finn Naes
% Date: April 14, 2026

clear;clc;close all

%% ASEN 3802 Lab 3 Part 1: Analysis of 2D Airfoils

%% Task 1 NACA 4-Digit Airfoil Generator

c=1;
N=50;

% NACA 0021
[x1,y1,xc_0021,yc_0021] = NACA_Airfoils(0,0,0.21,c,N); % NACA 0021

% NACA 2421
[x2,y2,xcam,ycam] = NACA_Airfoils(0.02,0.4,0.21,c,N); % NACA 2421

figure
hold on
grid on
plot(x1, y1, 'b.-', 'MarkerSize',8)
plot(x2, y2, 'r.-', 'MarkerSize',8)
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
grid on
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

N = 50; % 100 total, 50 per side

alpha = linspace(-15,15,100);
chord = 1;

cl_0006 = zeros(1,length(alpha));
cl_0012 = zeros(1,length(alpha));
cl_0018 = zeros(1,length(alpha));

% Vortex Panel Method
for i = 1:length(alpha)
    [x1,y1,~,~] = NACA_Airfoils(0,0,0.06,chord,N);
    cl_0006(i) = Vortex_Panel(x1,y1,alpha(i));

    [x2,y2,~,~] = NACA_Airfoils(0,0,0.12,chord,N);
    cl_0012(i) = Vortex_Panel(x2,y2,alpha(i));
    
    [x3,y3,~,~] = NACA_Airfoils(0,0,0.18,chord,N);
    cl_0018(i) = Vortex_Panel(x3,y3,alpha(i));
end

% Thin Airfoil Theory
cl_tat = 2*pi*(alpha*pi/180);

% Pulling Experimental Data
exp_0006 = readmatrix('Airfoil0006Data.csv');
exp_0012 = readmatrix('Airfoil0012Data.csv');
alpha_exp_0006 = exp_0006(:,1); cl_exp_0006 = exp_0006(:,2);
alpha_exp_0012 = exp_0012(:,1); cl_exp_0012 = exp_0012(:,2);

% plot
figure
hold on
grid on
plot(alpha,cl_0006,'b-')
plot(alpha,cl_0012,'r-')
plot(alpha,cl_0018,'g-')

plot(alpha,cl_tat,'k--')

plot(alpha_exp_0006,cl_exp_0006,'bo')
plot(alpha_exp_0012,cl_exp_0012,'ro')
legend('NACA 0006','NACA 0012','NACA 0018','Thin Airfoil Theory','NACA 0006 Exp','NACA 0012 Exp','Location','best')
xlabel('Angle of Attack \alpha [Deg]','Interpreter','tex')
ylabel('Sectional Lift COefficient C_l','Interpreter','tex')
title('Effect of Airfoil Thickness on Lift Coefficient')
ylim([-5 5])
hold off

% tables
lin_idx = (alpha >= -3) & (alpha <= 3);

fit_0006 = polyfit(alpha(lin_idx),cl_0006(lin_idx),1);
fit_0012 = polyfit(alpha(lin_idx),cl_0012(lin_idx),1);
fit_0018 = polyfit(alpha(lin_idx),cl_0018(lin_idx),1);

slope_0006 = fit_0006(1); a0_0006 = -fit_0006(2)/fit_0006(1);
slope_0012 = fit_0012(1); a0_0012 = -fit_0012(2)/fit_0012(1);
slope_0018 = fit_0018(1); a0_0018 = -fit_0018(2)/fit_0018(1);

slope_tat = 2*pi*(pi/180);
a0_tat = 0;

lin_exp6 = (alpha_exp_0006 >= -3) & (alpha_exp_0006 <= 3);
fit_exp6 = polyfit(alpha_exp_0006(lin_exp6),cl_exp_0006(lin_exp6),1);
slope_exp6 = fit_exp6(1); a0_exp6 = -fit_exp6(2)/fit_exp6(1); % 0006 experimental data slope

lin_exp12 = (alpha_exp_0012 >= -3) & (alpha_exp_0012 <= 3);
fit_exp12 = polyfit(alpha_exp_0012(lin_exp12),cl_exp_0012(lin_exp12),1);
slope_exp12 = fit_exp12(1); a0_exp12 = -fit_exp12(2)/fit_exp12(1); % 0012 experimental data slope

slope_exp18 = NaN; a0_exp18 = NaN; % no experimental data for 0018

airfoils = {'NACA 0006';'NACA 0012';'NACA 0018'};

% AoA
disp('Zero Lift Angle of Attack (Degrees)')
table_aoa = table(airfoils,[a0_0006;a0_0012;a0_0018], ...
    [a0_tat;a0_tat;a0_tat],[a0_exp6;a0_exp12;a0_exp18], ...
    'VariableNames',{'Airfoil','Vortex_Panel_deg','Thin_Airfoil_Theory_deg','Experimental_deg'});
disp(table_aoa)

% C_L
disp('Sectional Lift Coefficient (1/Degree)')
table_slope = table(airfoils,[slope_0006;slope_0012;slope_0018], ...
    [slope_tat;slope_tat;slope_tat],[slope_exp6;slope_exp12;slope_exp18], ...
    'VariableNames',{'Airfoil','Vortex_Panel','Thin_Airfoil_Theory','Experimental'});
disp(table_slope)

%% Task 4 Effect of Airfoil Camber on Lift

% pull experimental data
Airfoil2412Data=readmatrix("Airfoil2412Data.csv");
Airfoil4412Data=readmatrix("Airfoil4412Data.csv");


N = ceil(xid);

alpha_vals = linspace(-100,100,500);
chord = 1;
cl1 = zeros(1,length(alpha_vals));
cl2 = zeros(1,length(alpha_vals));
cl3 = zeros(1,length(alpha_vals));

% Run VPM
for i = 1:length(alpha_vals)

    [x1,y1,xc1,yc1] = NACA_Airfoils(0,0,0.12,chord,N);
    cl1(i) = Vortex_Panel(x1,y1,alpha_vals(i));

    [x2,y2,xc2,yc2] = NACA_Airfoils(0.02,0.4,0.12,chord,N);
    cl2(i) = Vortex_Panel(x2,y2,alpha_vals(i));
    
    [x3,y3,xc3,yc3] = NACA_Airfoils(0.04,0.4,0.12,chord,N);
    cl3(i) = Vortex_Panel(x3,y3,alpha_vals(i));

end

% Find zero-lift AOA
c_intersect = zeros(1,length(alpha_vals));
x_axis = 1:length(alpha_vals);

[xid1,yid1] = polyxpoly(x_axis,cl1,x_axis,c_intersect);
[xid2,yid2] = polyxpoly(x_axis,cl2,x_axis,c_intersect);
[xid3,yid3] = polyxpoly(x_axis,cl3,x_axis,c_intersect);

zero_lift_AOA1 = alpha_vals(floor(xid1));
zero_lift_AOA2 = alpha_vals(ceil(xid2));
zero_lift_AOA3 = alpha_vals(ceil(xid3));

alpha_zero1 = interp1(1:length(alpha_vals), alpha_vals, xid1);
alpha_zero2 = interp1(1:length(alpha_vals), alpha_vals, xid2);
alpha_zero3 = interp1(1:length(alpha_vals), alpha_vals, xid3);

cl_tat = 2*pi*(alpha_vals*pi/180);


% Plot all lift slopes
figure();
hold on
plot(alpha_vals,cl1, 'Color', 'b', 'LineWidth', 1)
plot(alpha_vals,cl2, 'Color', 'r', 'LineWidth', 1)
plot(alpha_vals,cl3, 'Color', 'magenta', 'LineWidth', 1)

plot(exp_0012(:,1),exp_0012(:,2),'b--', 'LineWidth', 2)
plot(Airfoil2412Data(:,1),Airfoil2412Data(:,2),'r--', 'LineWidth', 2)
plot(Airfoil4412Data(:,1),Airfoil4412Data(:,2),'--','Color', 'magenta', 'LineWidth', 2)

plot(alpha_vals,cl_tat,'k--','LineWidth',1.5)

xlabel('AOA \alpha [Deg]', 'Interpreter', 'tex')
ylabel('C_{L}')
legend('NACA 0012', 'NACA 2412', 'NACA 4412','NACA 0012 Experimental', 'NACA 2412 Experimental', 'NACA 4412 Experimental', 'Thin Airfoil Theory', 'Location', 'Best')
title('C_{L} vs. \alpha', 'Interpreter', 'tex')

xlim([-25,25])
ylim([-2,2])

hold off

zero_lift_angles = [alpha_zero1; alpha_zero2; alpha_zero3];
airfoil_names = {'NACA 0012'; 'NACA 2412'; 'NACA 4412'};


Vortex_Results = [alpha_zero1; alpha_zero2; alpha_zero3];
TAT_Results = [0; -2.08; -4.16];
Exp_Results = [0; -2.1; -4.3];

% Compare zero-lift AOA 
ComparisonTable = table(airfoil_names, Vortex_Results, TAT_Results, Exp_Results, 'VariableNames', {'Airfoil', 'Vortex_Panel_deg', 'Thin_Airfoil_Theory_deg', 'Experimental_deg'});
disp(ComparisonTable);

% Calculate lift slopes
linear_portion_VPM = (alpha_vals >= -3) & (alpha_vals <= 3);
linear_portion_exp1 = (exp_0012(:,1) >= -3) & (exp_0012(:,1) <= 3);
linear_portion_exp2 = (Airfoil2412Data(:,1) >= -3) & (Airfoil2412Data(:,1) <= 3);
linear_portion_exp3 = (Airfoil4412Data(:,1) >= -3) & (Airfoil4412Data(:,1) <= 3);

fit1_VPM = polyfit(alpha_vals(linear_portion_VPM), cl1(linear_portion_VPM), 1);
fit2_VPM = polyfit(alpha_vals(linear_portion_VPM), cl2(linear_portion_VPM), 1);
fit3_VPM = polyfit(alpha_vals(linear_portion_VPM), cl3(linear_portion_VPM), 1);

fit1_exp = polyfit(exp_0012(linear_portion_exp1,1), exp_0012(linear_portion_exp1,2),1);
fit2_exp = polyfit(Airfoil2412Data(linear_portion_exp2,1), Airfoil2412Data(linear_portion_exp2,2), 1);
fit3_exp = polyfit(Airfoil4412Data(linear_portion_exp3,1), Airfoil4412Data(linear_portion_exp3,2), 1);

lift_slope1_VPM = fit1_VPM(1);
lift_slope2_VPM = fit2_VPM(1);
lift_slope3_VPM = fit3_VPM(1);

slope_theoretical = (2 * pi) * pi / 180;

lift_slope1_exp = fit1_exp(1);
lift_slope2_exp = fit2_exp(1);
lift_slope3_exp = fit3_exp(1);

% Compare Lift Slopes
SlopeTable = table(airfoil_names, [lift_slope1_VPM; lift_slope2_VPM; lift_slope3_VPM], ...
                   ones(3,1)*slope_theoretical, [lift_slope1_exp; lift_slope2_exp; lift_slope3_exp], ...
                   'VariableNames', {'Airfoil', 'Vortex_Panel', 'Thin_Airfoil_Theory', 'Experimental'});

disp(SlopeTable);

%% ASEN 3802 Lab 3 Part 2: Analysis of Finite Wings with Thick Airfoils

b = 10; 
a0_t = 2*pi; 
a0_r = 2*pi; 
c_r = 1;
N = 50; 
aero_t = 3; 
aero_r = 3;
geo_t = 4; 
geo_r = 4; 

AR_list = [4 6 8 10];
taper_ratio = linspace(0,1,100);
delta_all = zeros(length(AR_list), length(taper_ratio));

for k = 1:length(AR_list)
    AR = AR_list(k);

    for i = 1:length(taper_ratio)
        taper = taper_ratio(i); 
        c_t = taper * c_r;
        b = AR * (c_r + c_t) / 2; % AR = b^2/S ++ S = (cr+ct)b/2
        [e, c_L, c_Di] = PLLT(b, a0_t, a0_r, c_t, c_r, aero_t, aero_r, geo_t, geo_r, N);
        delta_all(k,i) = (1/e) - 1;
    end
end


plot(taper_ratio, delta_all, 'LineWidth',2)
title('$\delta$ vs. $\frac{c_t}{c_r}$','Interpreter', 'latex')
xlabel('$\frac{c_t}{c_r}$','Interpreter', 'latex');
ylabel('$\delta$','Interpreter', 'latex');
legend(compose('AR = %g', AR_list), 'Location', 'best');

%% ASEN 3802 Lab 3 Part 3: Analysis of Approximate Cessna 140 Wing Performance

b = 33 + 4/12; % span
c_r = 5 + 4/12; % root chord
c_t = 3 + 8.5/12; % tip chord
S = b*(c_r+c_t)/2; % wing area

% Vortex Panel Method -> find a0 and alpha_L=0
N = N*2;
alpha_test = [-2,0,2];

% tip (NACA 0012)
cl_tip = zeros(1,3);
for i = 1:3
    [x,y,~,~]=NACA_Airfoils(0,0,0.12,c,N);
    cl_tip(i) = Vortex_Panel(x,y,alpha_test(i));
end
p_tip = polyfit(alpha_test,cl_tip,1);
a0_t_deg = p_tip(1);
a0_t = rad2deg(a0_t_deg); % convert to 1/rad for PLLT
aero_t = -p_tip(2)/p_tip(1); % zero-lift aoa in degrees

% root (NACA 2412)
cl_root = zeros(1,3);
for i = 1:3
    [x,y,~,~] = NACA_Airfoils(0.02,0.4,0.12,c,N);
    cl_root(i) = Vortex_Panel(x,y,alpha_test(i));
end
p_root = polyfit(alpha_test,cl_root,1);
a0_r_deg = p_root(1);
a0_r = rad2deg(a0_r_deg); % convert to 1/rad for PLLT
aero_r = -p_root(2)/p_root(1); % zero-lift aoa in degrees

% geometric twist
geo_r = 1;
geo_t = 0;

% For deliverables 1 + 2
alpha_wing = 4;
geo_r_conv = geo_r + alpha_wing;
geo_t_conv = geo_t + alpha_wing;

max_terms = 50;
N_terms = 1:max_terms;
CL_conv = zeros(1,max_terms);
CDi_conv = zeros(1,max_terms);

for i = 1:max_terms
    [~,CL_conv(i),CDi_conv(i)] = PLLT(b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,geo_t_conv,geo_r_conv,N_terms(i));
end

% assume N = max_terms is exact solution for error calculation
CL_exact = CL_conv(end);
CDi_exact = CDi_conv(end);

rel_err_CL = abs(CL_conv - CL_exact) / CL_exact * 100;
rel_err_CDi = abs(CDi_conv - CDi_exact) / CDi_exact * 100;

err_thresholds = [10, 1, 0.1];
terms_req_CL = zeros(1,3);
terms_req_CDi = zeros(1,3);

for k = 1:3
    terms_req_CL(k) = find(rel_err_CL <= err_thresholds(k),1,'first');
    terms_req_CDi(k) = find(rel_err_CDi <= err_thresholds(k),1,'first');
end

% DELIVERABLE 1
disp('Convergence Table')
Conv_Table = table(err_thresholds',terms_req_CL',terms_req_CDi',CL_conv(terms_req_CL)',CDi_conv(terms_req_CDi)','VariableNames',{'Error_Threshold_Percent','Odd_Terms_CL','Odd_Terms_CDi','CL_Value','CDi_Value'});
disp(Conv_Table)

% DELIVERABLE 2
figure

subplot(2,1,1)
plot(N_terms,CL_conv,'b-','LineWidth',2)
hold on
grid on
xline(terms_req_CL(1),'r--','10% Error')
xline(terms_req_CL(2),'g--','1% Error')
xline(terms_req_CL(3),'k--','0.1% Error')
xlabel('Number of Odd Terms'); ylabel('C_L','Interpreter','tex')
title('Convergence of Lift Coefficient')

subplot(2,1,2)
plot(N_terms,CDi_conv,'b-','LineWidth',2)
hold on
grid on
xline(terms_req_CDi(1),'r--','10% Error')
xline(terms_req_CDi(2),'g--','1% Error')
xline(terms_req_CDi(3),'k--','0.1% Error')
xlabel('Number of Odd Terms'); ylabel('C_{D,i}','Interpreter','tex')
title('Convergence of Induced Drag Coefficient')

%% Deliverable 3

N_del3 = max(terms_req_CL(3), terms_req_CDi(3));

h = 10000 * 0.3048; % ft to m
rho_imperial = 0.001756;% slug/ft^3
V_kts = 100;
V_fps = V_kts * 1.68781; % knots to ft/s
q = 0.5 * rho_imperial * V_fps^2; % lb/ft^2

alpha_wing = 4;
geo_r_conv = geo_r + alpha_wing;
geo_t_conv = geo_t + alpha_wing;

[~, CL_del3, CDi_del3] = PLLT(b, a0_t, a0_r, c_t, c_r, aero_t, aero_r, geo_t_conv, geo_r_conv, N_del3);

cd_0012 = 0.006; % From Abbot_TheoryOfWingSections_1959
cd_2412 = 0.0065; % Also from Abbot_TheoryOfWingSections_1959
cd_avg  = (cd_0012 + cd_2412) / 2;

CD_total = cd_avg + CDi_del3;

S_ft2 = S;

L = q * S_ft2 * CL_del3;
Di = q * S_ft2 * CDi_del3;
D  = q * S_ft2 * CD_total;
LD = L / D;
 
Del3_Table = table(L, Di, D, LD,'VariableNames', {'Lift_L_lbf','Induced_Drag_Di_lbf', 'Total_Drag_D_lbf','L_over_D'});
disp('Deliverable 3 Results')
disp(Del3_Table)

%% Deliverable 4

alpha_vals = linspace(-5, 15, 100);
N_del4 = max(terms_req_CL(3), terms_req_CDi(3));

CL_sweep  = zeros(1, length(alpha_vals));
CDi_sweep = zeros(1, length(alpha_vals));

for i = 1:length(alpha_vals)
    geo_r_i = geo_r + alpha_vals(i);
    geo_t_i = geo_t + alpha_vals(i);

    [~, CL_sweep(i), CDi_sweep(i)] = PLLT(b, a0_t, a0_r, c_t, c_r, aero_t, aero_r, geo_t_i, geo_r_i, N_del4);
end

% [Cl, Cd] 
Abbott_0012 = [-0.8, 0.0080; -0.4, 0.0065; 0.0, 0.0060; 0.4, 0.0065; 0.8, 0.0080; 1.0, 0.0100; 1.2, 0.0130];
Abbott_2412 = [-0.4, 0.008; 0.0, 0.0065; 0.4, 0.006; 0.8, 0.0070; 1.0, 0.0085; 1.2, 0.0110; 1.4, 0.0140];

cd_0012_interp = interp1(Abbott_0012(:,1), Abbott_0012(:,2), CL_sweep, 'linear', 'extrap');
cd_2412_interp = interp1(Abbott_2412(:,1), Abbott_2412(:,2), CL_sweep, 'linear', 'extrap');

cd_profile_interp = (cd_0012_interp + cd_2412_interp) / 2;

CD_total_interp = cd_profile_interp + CDi_sweep;

figure();
hold on
grid on

% Plot all components
plot(alpha_vals, CD_total_interp, 'k-', 'LineWidth', 2)
plot(alpha_vals, CDi_sweep, 'r--', 'LineWidth', 2)
plot(alpha_vals, cd_profile_interp,'b--', 'LineWidth', 2)
xlabel('Angle of Attack \alpha [deg]', 'Interpreter', 'tex')
ylabel('Drag Coefficient C_D', 'Interpreter', 'tex')
title('Total Drag Coefficient vs. Angle of Attack')
legend('C_D Total', 'C_{D,i} Induced', 'c_d Profile', 'Location', 'best', 'Interpreter', 'tex')
hold off

%% Deliverable 5

L_D_total = CL_sweep./CD_total_interp;

figure();
hold on
grid on
plot(alpha_vals,L_D_total,'k','LineWidth',2)

xlabel('Angle of Attack \alpha [deg]')
ylabel('Lift to Drag Ratio L/D')
title('Lift to Drag Ratio vs. Angle of Attack')

hold off




%% FUNCTIONS

function [x_b, y_b, x, y_c] = NACA_Airfoils(m,p,t,c,N)
% Constructs panels for any NACA 4-digit airfoil.
% 
% Inputs: max camber (m), max camber location (p), max thickness (t), chord
% length (c), number of panels (N).
% 
% Calculates shape by determining mean camber and applying thickness
% distribution perpendicular to camber.
%
% Clusters points more tightly at LE and TE using equiangular spacing.
%
% Outputs: x_b, y_b (x- and y- coords of airfoils surface boundary points)

    % cosine (equiangular spacing)
    beta = linspace(0,pi,N+1);
    x = (c/2)*(1-cos(beta));

    y_t = (t/0.2)*c*(0.2969*sqrt(x/c)-0.1260*(x/c)-0.3516*(x/c).^2+0.2843*(x/c).^3-0.1036*(x/c).^4);

    y_c = zeros(size(x));
    dyc_dx = zeros(size(x));

    for i = 1:length(x)
        if x(i) < p*c
            if p == 0
                y_c(i) = 0;
                dyc_dx(i) = 0;
            else
                y_c(i) = m*(x(i)/p^2)*(2*p-x(i)/c);
                dyc_dx(i) = (2*m/p^2)*(p-x(i)/c);
            end
        else
            if p == 1
                y_c(i) = 0;
                dyc_dx(i) = 0;
            else
                y_c(i) = m*(c-x(i))/(1-p)^2*(1+x(i)/c-2*p);
                dyc_dx(i) = (2*m/(1-p)^2)*(p-x(i)/c);
            end
        end
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

function [e,c_L,c_Di] = PLLT(b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,geo_t,geo_r,N)
    % Description: Solves the fundamental equation of PLLT for finite wings
    % with thick airfoils.
    %
    % Inputs:
    % span (b)
    % tip cross-sectional lift slope (a0_t)
    % root cross-sectional lift slope (a0_r)
    % tip chord (c_t), root chord (c_r)
    % tip zero-lift aoa (aero_t), root zero-lift aoa (aero_r)
    % tip geometric aoa (geo_t), root geometric aoa (geo_r)
    % odd terms included in series expansion for circulation (N)
    %
    % Outputs:
    % span efficiency factor (e)
    % coefficient of lift (c_L)
    % induced coefficient of drag (c_Di)
    % degrees to radians

    aero_t = deg2rad(aero_t);
    aero_r = deg2rad(aero_r);
    geo_t = deg2rad(geo_t);
    geo_r = deg2rad(geo_r);

    % initialize M and RHS
    M = zeros(N,N);
    RHS = zeros(N,1);

    % build system of equations at N locations
    for i = 1:N
        theta_i = (i*pi)/(2*N); % eq 2
        eta = abs(cos(theta_i)); % spanwise interpolation factor (eta)

        % linearly interpolate properties at station i
        a0_i = a0_r + (a0_t - a0_r)*eta;
        c_i = c_r + (c_t - c_r)*eta;
        aero_i = aero_r + (aero_t - aero_r)*eta;
        geo_i = geo_r + (geo_t - geo_r)*eta;

        RHS(i) = geo_i - aero_i; % RHS of PLLT equation

        % eq 1, 3
        for j = 1:N
            n = 2*j - 1; % odd terms

            term1 = (4*b)/(a0_i*c_i);
            term2 = n/sin(theta_i);

            M(i,j) = sin(n*theta_i)*(term1+term2); % M matrix
        end
    end

    A = M\RHS; % linear system for fourier coefficients (A_n)

    S = ((c_r + c_t)/2)*b; % wing area

    AR = (b^2)/S;% aspect ratio

    A1 = A(1);
    c_L = A1*pi*AR; % lift coefficient

    c_Di = 0; % induced drag coefficient
    for j = 1:N
        n = 2*j-1;
        c_Di = c_Di + pi*AR*n*(A(j)^2);
    end

    e = (c_L^2)/(pi*AR*c_Di); % span efficiency factor
end
