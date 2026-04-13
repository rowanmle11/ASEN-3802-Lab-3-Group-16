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