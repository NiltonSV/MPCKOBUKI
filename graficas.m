close all
clear
clc

data = readmatrix('results.txt'); % Cargar datos sin problemas con encabezados

x0 = data(:,1);
x1 = data(:,2);
x2 = data(:,3);
delta_qr0 = data(:,4);
delta_qr1 = data(:,5);

figure;
subplot(3,1,1);
hold on;
plot(x0, x1, 'b-');
% plot(1:length(x1), x2)
title('Trayectoria del Robot Real');
grid on;

subplot(3,1,2);
hold on;
plot(1:length(x1), x2)
title('Orientación del Robot Real');
grid on;

subplot(3,1,3);
plot(delta_qr0, 'r-', 'DisplayName', 'v');
hold on;
plot(delta_qr1, 'g-', 'DisplayName', 'w');
xlabel('Iteraciones');
ylabel('Valores de velocidades Reales');
title('Entradas de Control Reales');
legend;
grid on;

data = readmatrix('reference.txt'); % Cargar datos sin problemas con encabezados

xr0 = data(:,1);
xr1 = data(:,2);
xr2 = data(:,3);
deltar_qr0 = data(:,4);
deltar_qr1 = data(:,5);

% hold on
% th = 0:pi/50:2*pi;
% xunit = r * cos(th) + x;
% yunit = r * sin(th) + y;

figure;
subplot(3,1,1);
plot(xr0, xr1, 'b-', 'LineWidth', 2);
hold on;
plot(x0, x1, 'r-', 'LineWidth', 2);
hold on;
dataobst = readmatrix('obstacles.txt');
if ~all(isnan(dataobst), 'all')
    for i=1:size(dataobst,1)
        thc = 0:pi/50:2*pi;
        xc = dataobst(i,3) * cos(thc) + dataobst(i,1);
        yc = dataobst(i,3) * sin(thc) + dataobst(i,2);
        plot(xc, yc, 'black-', 'LineWidth', 2);
        hold on;
        xc = (dataobst(i,3)+0.02) * cos(thc) + dataobst(i,1);
        yc = (dataobst(i,3)+0.02) * sin(thc) + dataobst(i,2);
        plot(xc, yc, 'y-', 'LineWidth', 2);
        hold on;
    end
end
% thc = 0:pi/50:2*pi;
% xc = 0.3 * cos(thc) + 1;
% yc = 0.3 * sin(thc) - 0;
% plot(xc, yc, 'y-');
% hold on;
% plot(1:length(xr0), xr0 - x0)
title('xye');
% xlim([-1.5 1.5]);
% ylim([-2 2]);
grid on;

subplot(3,1,2);
plot(1:length(xr0),xr0 - x0, 'b-');
hold on;
plot(1:length(xr1), xr1 - x1)
title('ye xe');
legend('xe', 'ye')
grid on;

% Valores de velocidades
subplot(3,1,3);
plot(1:length(xr2),xr2 - x2, 'r-');
hold on;
% plot(deltar_qr0, 'r-', 'DisplayName', 'v');
% hold on;
% plot(deltar_qr1, 'g-', 'DisplayName', 'w');
% xlabel('Iteraciones');
title('Th e');
% title('Entradas de Control');
legend('th');
grid on;

data = readmatrix('reference.txt'); % Cargar datos sin problemas con encabezados

x0 = data(:,1);
x1 = data(:,2);
x2 = data(:,3);
delta_qr0 = data(:,4);
delta_qr1 = data(:,5);

figure;
subplot(3,1,1);
hold on;
plot(x0, x1, 'b-');
% plot(1:length(x1), x2)
title('Trayectoria del Robot de Referencia');
grid on;

subplot(3,1,2);
hold on;
plot(1:length(x1), x2)
title('Orientación del Robot de Referencia');
grid on;

subplot(3,1,3);
plot(delta_qr0, 'r-', 'DisplayName', 'vref');
hold on;
plot(delta_qr1, 'g-', 'DisplayName', 'wref');
xlabel('Iteraciones');
ylabel('Valores de velocidades');
title('Entradas de Control');
legend;
grid on;
