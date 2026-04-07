function [x_b, y_b, x, y_c] = NACA_Airfoils(m,p,t,c,N)

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
    yL = -(y_c+y_t.*cos(xi));

    x_lower = xL(end:-1:1);
    y_lower = yL(end:-1:1);

    x_upper = xU(2:end);
    y_upper = yU(2:end);

    x_b = [x_lower,x_upper];
    y_b = [y_lower,y_upper];
end