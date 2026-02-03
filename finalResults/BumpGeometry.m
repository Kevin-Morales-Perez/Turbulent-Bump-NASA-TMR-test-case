%Geometry of the bump
%{
x=linspace(0,1.5,100);
%y=0.05*sin(pi*x/0.9 -(pi/3)).^4;
y=zeros(1,100);

for i=1:100
    if (x(i)<=0.3) || (x(i)>=1.2)
        y(i)=0;
    else
        y(i)=0.05*sin(pi*x(i)/0.9 -(pi/3)).^4;
        
    end
end

plot(x,y)
title('Detalle de protuberancia para caso de prueba presentado en el \itTurbulence Modeling Resource \rm','LineWidth',6)
xlabel('X (m)')
ylabel('Y (m)')
grid on
%}

% Geometry of the bump (NASA TMR test case)

x = linspace(0,1.5,400);
y = zeros(size(x));

for i = 1:length(x)
    if (x(i) <= 0.3) || (x(i) >= 1.2)
        y(i) = 0;
    else
        y(i) = 0.05 * sin(pi*x(i)/0.9 - pi/3).^4;
    end
end

figure
plot(x, y, 'b', 'LineWidth', 2.5)
title('Detalle de protuberancia para caso de prueba presentado en el \itTurbulence Modeling Resource \rm ', ...
       'FontSize', 12)

xlabel('x (m)')
ylabel('y (m)')
grid on
xlim([0 1.5])
ylim([0 0.06])
