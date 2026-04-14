clc;
clear;
close all;
b = 10; 
a0_t = 2*pi; 
a0_r = 2*pi; 
c_r = 1;
N = 50; 
aero_t = 3; 
aero_r = 3;
geo_t = 4; 
geo_r = 4; 

AR_Vals = [4 6 8 10];
taper_ratio = linspace(0.008,1,100);
delta_all = zeros(length(AR_Vals), length(taper_ratio));

for k = 1:length(AR_Vals)
    AR = AR_Vals(k);
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
legend(compose('AR = %g', AR_Vals), 'Location', 'best');

function [e,c_L,c_Di] = PLLT(b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,geo_t,geo_r,N)

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
