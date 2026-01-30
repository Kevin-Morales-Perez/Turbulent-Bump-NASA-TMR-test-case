%% BLOCK 1 - extracting data from dat file
clear; clc; close all;


filename = 'mut_contours_cfl3d.dat';

% Open and read entire file
fid = fopen(filename, 'r');

% Read header lines manually
for i = 1:8
    fgetl(fid);
end

% The 9th line contains the data type info
line9 = fgetl(fid);
if ~contains(line9, 'DT=')
    error('Unexpected file format');
end

% Now read all numbers from the rest of the file
all_numbers = textscan(fid, '%f', 'CollectOutput', true);
fclose(fid);

all_numbers = all_numbers{1};  % Extract from cell array

% Parameters
I = 1408;
J = 640;
total_pts = I * J;

% Split
X_data = all_numbers(1:total_pts);
Y_data = all_numbers(total_pts+1:2*total_pts);
mut_data = all_numbers(2*total_pts+1:3*total_pts);

% Reshape
X_grid = reshape(X_data, [I, J])';
Y_grid = reshape(Y_data, [I, J])';
mut_grid = reshape(mut_data, [I, J])';

%% plottin own code

%load("finalConvergedRe3e6.mat")
load("ReSolveBoundTurbItNum_900000.mat")
%levels_nondimEddyvisc=linspace(0,600,10);%use for specific levels

%for i=ny:-1:1
%y_temp=Yctrs(i,1);
%if y_temp>0.08
%disp(i)
%break
%end
%end
%    15

mu_nondim=mu_turbulent/mu;
x_zoom_bump=Xctrs(:,nxSolid);
y_zoom_bump=Yctrs(:,nxSolid);
mu_nondim_zoom_bump=mu_nondim(:,nxSolid);
levels_nondimEddyvisc=linspace(0,260,21);%% correct levels !

figure(1)%Non Dimensional Eddy viscosity
contourf(x_zoom_bump,y_zoom_bump,mu_nondim_zoom_bump,levels_nondimEddyvisc , "ShowText","on")
title("Non dimensional Eddy Viscosity ")
xlabel("Lenght (m)")
ylabel("Height (m)")
colormap jet
colorbar
axis([distPlat distPlat+1.5 0 0.08])
%{

figure(1)%Non dimensional Eddy viscosity  own code
contourf(Xctrs(nxSolid,ny:-1,10),Yctrs(nxSolid,ny:-1,10),mu_turbulent/mu, levels_nondimEddyvisc, "ShowText","on")
title("Non dimensional eddy viscosity")
xlabel("Lenght (m)")
ylabel("Height (m)")
colormap jet
colorbar
axis([distPlat distPlat+1.5 0 0.08])
%}

%% Ploting NASA cfl3 results 
figure(2)%Non Dimensional Eddy viscosity
contourf(X_grid,Y_grid,mut_grid,  levels_nondimEddyvisc, "ShowText","on")
title("CFL3 Non dimensional Eddy Viscosity ")
xlabel("Lenght (m)")
ylabel("Height (m)")
colormap jet
colorbar
axis([0 1.5 0 0.08])

%% comparing both results 


%% Side-by-side with shared colorbar

% Left plot
figure(3)
subplot(1, 2, 1);
contourf(x_zoom_bump,y_zoom_bump,mu_nondim_zoom_bump,levels_nondimEddyvisc , "ShowText","on")
axis([distPlat distPlat+1.5 0 0.08])
% Right plot
subplot(1, 2, 2);
contourf(X_grid,Y_grid,mut_grid,  levels_nondimEddyvisc, "ShowText","on")
colormap jet
colorbar
axis([0 1.5 0 0.08])
%}

%% other plot

% Simplified corrected version
% Clean version without isoline labels (professional for publications)
figure(3)
set(gcf, 'Position', [100, 100, 1200, 500]);

% Left plot (Your simulation)
subplot(1, 2, 1);
contourf(x_zoom_bump, y_zoom_bump, mu_nondim_zoom_bump, levels_nondimEddyvisc, 'LineWidth', 0.5);
xlabel('Longitud (m)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Altura (m)', 'FontSize', 11, 'FontWeight', 'bold');
title('Simulación CFD Actual', 'FontSize', 12, 'FontWeight', 'bold');
axis([distPlat distPlat+1.5 0 0.08]);
grid on;

% Add Greek letter label
text(0.02, 0.85, '$\frac{\mu_t}{\mu}$', 'Units', 'normalized', ...
     'FontSize', 14, 'Interpreter', 'latex', 'Color', 'white', ...
     'BackgroundColor', [0.2 0.2 0.2 0.5], 'EdgeColor', 'white');

% Right plot (Reference data)
subplot(1, 2, 2);
contourf(X_grid, Y_grid, mut_grid, levels_nondimEddyvisc, 'LineWidth', 0.5);
xlabel('Longitud (m)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Altura (m)', 'FontSize', 11, 'FontWeight', 'bold');
title('Datos de Referencia', 'FontSize', 12, 'FontWeight', 'bold');
axis([0 1.5 0 0.08]);
grid on;

% Add Greek letter label
text(0.02, 0.85, '$\frac{\mu_t}{\mu}$', 'Units', 'normalized', ...
     'FontSize', 14, 'Interpreter', 'latex', 'Color', 'white', ...
     'BackgroundColor', [0.2 0.2 0.2 0.5], 'EdgeColor', 'white');

% Common colorbar
colormap(jet);
h = colorbar('Position', [0.93 0.15 0.02 0.7]);
h.Label.String = '$\frac{\mu_t}{\mu}$';
h.Label.Interpreter = 'latex';
h.Label.FontSize = 14;
h.Label.FontWeight = 'bold';

% Overall title
sgtitle('Viscosidad de Torbellino Adimensional $\left(\frac{\mu_t}{\mu}\right)$ a $Re = 3 \times 10^6$', ...
        'FontSize', 16, 'FontWeight', 'bold', 'Interpreter', 'latex');

%% 2 plot

% Enhanced clean version for thesis/publication
figure(3)
set(gcf, 'Position', [100, 100, 1300, 550], 'Color', 'white');

% Use a better colormap (optional - 'parula' is perceptually uniform)
% colormap(parula); % Uncomment for a more professional colormap

% Left plot
ax1 = subplot(1, 2, 1);
hold on;
contourf(x_zoom_bump, y_zoom_bump, mu_nondim_zoom_bump, levels_nondimEddyvisc, ...
         'LineColor', [0.3 0.3 0.3], 'LineWidth', 0.5);

% Add the bump geometry for context (if you have it)
% plot(x_bump, y_bump, 'k-', 'LineWidth', 2);

xlabel('Longitud (m)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Altura (m)', 'FontSize', 12, 'FontWeight', 'bold');
title('Simulación CFD Actual', 'FontSize', 13, 'FontWeight', 'bold');
axis([distPlat distPlat+1.5 0 0.08]);
grid on;
box on;
set(gca, 'FontSize', 11, 'GridAlpha', 0.3);

% Add Greek letter annotation in a box
annotation('textbox', [0.14, 0.77, 0.08, 0.08], ...
           'String', '$\frac{\mu_t}{\mu}$', ...
           'Interpreter', 'latex', ...
           'FontSize', 16, ...
           'FontWeight', 'bold', ...
           'EdgeColor', 'black', ...
           'BackgroundColor', 'white', ...
           'FitBoxToText', 'on');

% Right plot
ax2 = subplot(1, 2, 2);
hold on;
contourf(X_grid, Y_grid, mut_grid, levels_nondimEddyvisc, ...
         'LineColor', [0.3 0.3 0.3], 'LineWidth', 0.5);

xlabel('Longitud (m)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Altura (m)', 'FontSize', 12, 'FontWeight', 'bold');
title('Datos de Referencia (NASA TMR)', 'FontSize', 13, 'FontWeight', 'bold');
axis([0 1.5 0 0.08]);
grid on;
box on;
set(gca, 'FontSize', 11, 'GridAlpha', 0.3);

% Add Greek letter annotation in a box
annotation('textbox', [0.64, 0.77, 0.08, 0.08], ...
           'String', '$\frac{\mu_t}{\mu}$', ...
           'Interpreter', 'latex', ...
           'FontSize', 16, ...
           'FontWeight', 'bold', ...
           'EdgeColor', 'black', ...
           'BackgroundColor', 'white', ...
           'FitBoxToText', 'on');

% Common colorbar with better formatting
h = colorbar('Position', [0.92 0.15 0.018 0.7]);
h.Label.String = '$\frac{\mu_t}{\mu}$';
h.Label.Interpreter = 'latex';
h.Label.FontSize = 15;
h.Label.FontWeight = 'bold';
h.Label.Rotation = 0; % Horizontal label
h.Label.Position = [2.8, 0.5, 0]; % Adjust position
h.TickDirection = 'out';
h.LineWidth = 1.2;
h.FontSize = 11;

% Format tick labels with consistent decimal places
ticks = get(h, 'Ticks');
if ~isempty(ticks)
    tick_labels = arrayfun(@(x) sprintf('%.2f', x), ticks, 'UniformOutput', false);
    set(h, 'TickLabels', tick_labels);
end

% Main title
sgtitle('Comparación de Viscosidad de Torbellino Adimensional $\left(\frac{\mu_t}{\mu}\right)$ a $Re = 3 \times 10^6$', ...
        'FontSize', 17, 'FontWeight', 'bold', 'Interpreter', 'latex');

% Optional: Add a note about the bump location
annotation('textbox', [0.15, 0.05, 0.7, 0.04], ...
           'String', 'Nota: La protuberancia se extiende desde x = 0.5 m hasta x = 1.0 m aproximadamente', ...
           'FontSize', 10, 'HorizontalAlignment', 'center', ...
           'EdgeColor', 'none', 'BackgroundColor', 'none');

% Ensure both plots have the same colormap limits
caxis(ax1, [min(levels_nondimEddyvisc), max(levels_nondimEddyvisc)]);
caxis(ax2, [min(levels_nondimEddyvisc), max(levels_nondimEddyvisc)]);

%% 3 plot

% Minimalist clean version - no annotations, just the plots
figure(3)
set(gcf, 'Position', [100, 100, 1100, 450]);

% Left plot
subplot(1, 2, 1);
contourf(x_zoom_bump, y_zoom_bump, mu_nondim_zoom_bump, levels_nondimEddyvisc);
xlabel('Longitud (m)');
ylabel('Altura (m)');
title('Simulación CFD');
axis([distPlat distPlat+1.5 0 0.08]);

% Right plot
subplot(1, 2, 2);
contourf(X_grid, Y_grid, mut_grid, levels_nondimEddyvisc);
xlabel('Longitud (m)');
ylabel('Altura (m)');
title('Referencia (NASA)');
axis([0 1.5 0 0.08]);

% Colorbar
colormap(jet);
colorbar;
grid on
% Title
sgtitle('Viscosidad de torbellino adimensional (\mu_t/\mu), Re=3\times10^6');

%% other plot 
% Minimalist clean version - no annotations, just the plots
figure(3)
set(gcf, 'Position', [100, 100, 1100, 450]);

% Left plot
subplot(1, 2, 1);
contourf(x_zoom_bump, y_zoom_bump, mu_nondim_zoom_bump, levels_nondimEddyvisc, ...
         'LineColor', 'k');  % Black isolines
xlabel('Longitud (m)');
ylabel('Altura (m)');
title("Código propio");
axis([distPlat distPlat+1.5 0 0.08]);
pbaspect([1 1 1]); 
grid on;

% Right plot
subplot(1, 2, 2);
contourf(X_grid, Y_grid, mut_grid, levels_nondimEddyvisc, ...
         'LineColor', 'k');  % Black isolines
xlabel('Longitud (m)');
ylabel('Altura (m)');
title('NASA CFL3');
axis([0 1.5 0 0.08]);
pbaspect([1 1 1]); 
grid on;

% Colorbar
colormap(jet);
colorbar;

% Title with Greek letters using LaTeX interpreter
sgtitle("Comparaci\'on de viscosidad de torbellino adimensional ($\mu_t/\mu$), $Re=3\times10^6$", ...
        'Interpreter', 'latex');

%% MY plot
figure(4)
subplot(1, 2, 1);%OWN RESULTS
contourf(x_zoom_bump, y_zoom_bump, mu_nondim_zoom_bump, levels_nondimEddyvisc, ...
         'LineColor', 'k')
xlabel('Longitud (m)')
ylabel('Altura (m)')
title("Código propio")
grid on
axis([distPlat distPlat+1.5 0 0.08])


subplot(1, 2, 2);%NASA RESULTS
contourf(X_grid, Y_grid, mut_grid, levels_nondimEddyvisc,'LineColor', 'k')  % Black isolines
xlabel('Longitud (m)')
ylabel('Altura (m)')
title('NASA CFL3D')
grid on
axis([0 1.5 0 0.08])

% Common colorbar
colormap(jet);
h = colorbar('Position', [0.93 0.15 0.02 0.7]);
h.Label.String = '$\frac{\mu_t}{\mu}$';
h.Label.Interpreter = 'latex';
h.Label.FontSize = 14;
h.Label.FontWeight = 'bold';

sgtitle("Comparaci\'on de viscosidad de torbellino adimensional ($\mu_t/\mu$), $Re=3\times10^6$", ...
        'Interpreter', 'latex')
