function [x_b,y_b] = NACA(m,p,t,c,N)

x = linspace(0,c,N+1);


y_t = (t/0.2)*c.*(0.2969.*sqrt(x./c)-0.1260.*(x./c)-0.3516.*(x./c).^2+0.2843.*(x./c).^3-0.1036.*(x./c).^4);



yc = zeros(1,length(x));
dy_dx = zeros(1,length(x));

for i = 1:(N+1)
    if x(i) < p * c
        yc(i) = m .* (x(i)./p^2).*((2*p)-x(i)./c);
        dy_dx(i) = (m / p^2) .* (2*p - x(i)./c)+(m.*x(i)./(p^2)).*(-1/c);
    else
        yc(i) = m .* ((c - x(i)) ./ (1-p)^2) .* (1 + x(i)./c - 2*p);
        dy_dx(i) = (-m / (1-p)^2)*(1+x(i)./c-2*p)+((m.*(c-x(i))./(1-p)^2).*(1/c));
    end
end

A = atan(dy_dx);

xu = x - y_t.*(sin(A));
xl = x + y_t.*(sin(A));

yu = yc + y_t.*(cos(A));
yl = yc - y_t.*(cos(A));



figure;
hold on
plot(x, y_t, "Color",'r');
plot(x,yc,"Color",'b');
plot(xu,yu,"Color",'k');
plot(xl,yl,"Color",'k');
ylim([-0.5,0.5])
hold off;
end
