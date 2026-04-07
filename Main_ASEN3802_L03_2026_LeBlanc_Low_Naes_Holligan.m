% Contributors: Rowan LeBlanc

clear;clc;close all

%% Task 1

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

%% Task 2



%% Task 3



%% Task 4