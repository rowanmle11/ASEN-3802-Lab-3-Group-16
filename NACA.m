function [x_b,y_b] = NACA(m,p,t,c,N)

x_angle = linspace(0,c,N);


y_t = (t/0.2)*c.*(0.2969.*sqrt(x_angle./c)-0.1260.*(x_angle./c)-0.3516.*(x_angle./c).^2+0.2843.*(x_angle./c).^3-0.1036.*(x_angle./c).^4);



yc = zeros(1,length(x_angle));
dy_dx = zeros(1,length(x_angle));

for i = 1:N
    if x_angle(i) < p * c
        yc(i) = m .* (x_angle(i)./p^2).*((2*p)-x_angle(i)./c);
        dy_dx(i) = (m / p^2) .* (2*p - x_angle(i)./c)+(m.*x_angle(i)./(p^2)).*(-1/c);
    else
        yc(i) = m .* ((c - x_angle(i)) ./ (1-p)^2) .* (1 + x_angle(i)./c - 2*p);
        dy_dx(i) = (-m / (1-p)^2)*(1+x_angle(i)./c-2*p)+((m.*(c-x_angle(i))./(1-p)^2).*(1/c));
    end
end

A = atan(dy_dx);

xu = x_angle - y_t.*(sin(A));
xl = x_angle + y_t.*(sin(A));

yu = yc + y_t.*(cos(A));
yl = yc - y_t.*(cos(A));



figure;
hold on
plot(x_angle, y_t, "Color",'r');
plot(x_angle,yc,"Color",'b');
plot(xu,yu,"Color",'k');
plot(xl,yl,"Color",'k');
ylim([-0.5,0.5])
hold off;
end