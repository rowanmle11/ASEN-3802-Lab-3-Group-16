function [x_b, y_b, x_camber, y_camber] = NACA_Airfoils(m,p,t,c,N)

    N_half = N/2;

    % cosine (equiangular spacing)
    theta = linspace(0,pi,N_half+1);
    x = (c/2)*(1-cos(theta));

    y_t = (t/0.2)*c*(0.2969*sqrt(x/c)-0.1260*(x/c)-0.3516*(x/c).^2+0.2843*(x/c).^3-0.1036*(x/c).^4);

    y_c = zeros(size(x));
    dyc_dx = zeros(size(x));

    for i = 1:length(x)
        if x(i) >= 0 && x(i) < p*c
            y_c = m.*(x./p.^2).*((2.*p)-(x./c));
            dyc_dx = -(2.*m.*(x-(c.*p)))./(c.*p.^2);
        elseif x(i) >= p*c && x(i) <= c
            y_c = m.*((c-x)./(1-p).^2).*(1+(x./c)-(2.*p));
            dyc_dx = -(2.*m.*(x-(c.*p)))./(c.*(p-1).^2);
        end
    end

    xi = atan(dyc_dx);

    % upper
    xU = x-y_t.*sin(xi);
    yU = y_c+y_t.*cos(xi);
    % lower
    xL = x+y_t.*sin(xi);
    yL = -(y_c+y_t.*cos(xi));

    x_b = [flip(xU),xL(2:end)];
    y_b = [flip(yU), yL(2:end)];

    x_camber = x;
    y_camber = y_c;
end