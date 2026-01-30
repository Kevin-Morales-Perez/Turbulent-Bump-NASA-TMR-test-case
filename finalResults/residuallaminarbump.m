step=15000; 

figure(2)
plot(1:step:iterations_cont-1,...
    residualsMat(1:step:iterations_cont-1,1),'-s',...
    1:step:iterations_cont-1,...
    residualsMat(1:step:iterations_cont-1,2),"-o",...
    1:step:iterations_cont-1,...
    residualsMat(1:step:iterations_cont-1,3),"-d",...
    1:step:iterations_cont-1,...
    residualsMat(1:step:iterations_cont-1,4),"-p",...
    1:step:iterations_cont-1,...
    residualsMat(1:step:iterations_cont-1,5),"->")

legend('Velocidad u','Velocidad v','Corrección de presión','Turbulencia','Continuidad')
title('Historial de convergencia de los residuales')
xlabel('Iteraciones')
ylabel('Residual (escala logarítmica)')
yscale log
grid on
pbaspect([2 1.25 1.25]) 
%% u pplot

figure(4)
contourf(Xctrs,Yctrs,u, 20, 'LineColor', 'none')
title('Componente u de la velocidad (m/s)')
xlabel('Longitud (m)')
ylabel('Altura (m)')
colormap jet
colorbar
axis equal
xlim([3 4.5])
ylim([0 1])

%% V
figure(5)
contourf(Xctrs,Yctrs,v, 20, 'LineColor', 'none') 
title('Componente v de la velocidad (m/s)') 
xlabel('Longitud (m)') 
ylabel('Altura (m)') 
colormap jet
colorbar
axis equal
xlim([2.75 4.75])
ylim([0 1])


%% P

figure(6)
contourf(Xctrs, Yctrs, p, 20, 'LineColor', 'none')
hold on
[C,h] = contour(Xctrs, Yctrs, p, 20, 'k', 'LineWidth', 0.6);
clabel(C,h,'FontSize',8,'Color','k')
hold off

title('Presión (Pa)')
xlabel('Longitud (m)')
ylabel('Altura (m)')
colormap jet
colorbar
axis equal
xlim([2.5 5])
ylim([0 1.5])
%% Variable transportada nu_tilde (Spalart–Allmaras)

figure(13)
contourf(Xctrs, Yctrs, nu_tilde, 20, 'LineColor', 'none')
title('Variable transportada $\tilde{\nu}$ del modelo de Spalart--Allmaras (m$^2$/s)', ...
      'Interpreter','latex')
xlabel('Longitud (m)')
ylabel('Altura (m)')
colormap(jet)
colorbar
axis equal
xlim([3 5.5])
ylim([0 1])



%% Viscosidad turbulenta mu_t

figure(14)
contourf(Xctrs, Yctrs, mu_turbulent, 20, 'LineColor', 'none')
title('Viscosidad dinámica turbulenta \mu_t (Pa·s)')
xlabel('Longitud (m)')
ylabel('Altura (m)')
colormap(jet)
colorbar
axis equal
xlim([3 5.5])
ylim([0 1])



