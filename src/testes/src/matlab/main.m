clc, clear, close all;

x = [0; 0.1; 0];
C = eye(3);
dt = 0.1; % sample time

% Path
traj = "u";
[xref, yref, thref, vref, wref] = Path(dt, traj);
n_points                        = length(xref);

% N and M
N = 15;
M = 3;

% Q and R matrices
Q = [1.0 0 0;
     0 1.0 0;
     0 0 0.5];

R = [0.1 0;
     0 0.1];

% Boundaries
lb = [-1.47; -3.77];
ub = [1.47; 3.77];
for i=2:M
    lb = [lb, [-1.47; -3.77]];
    ub = [ub, [1.47; 3.77]];
end

xrobot = zeros(1, n_points);
yrobot = zeros(1, n_points);
throbot = zeros(1, n_points);
vrobot  = zeros(1, n_points);
wrobot  = zeros(1, n_points);

for i=1:n_points

    xrobot(i) = x(1);
    yrobot(i) = x(2);
    throbot(i) = x(3);

    A = [1 0 -vref(i)*sin(thref(i))*dt;
        0 1 vref(i)*cos(thref(i))*dt;
        0 0 1];

    B = [cos(thref(i))*dt 0;
        sin(thref(i))*dt 0;
        0 dt];

    % **************************** Referencia ****************************

    % Ref = [];

    % for j=1:N
    %     if i + j - 1 < n_points
    Ref = [xref(i); yref(i); thref(i)];
    % Ref = [Ref; ref];
    %     else
    %         ref = [xref(n_points); yref(n_points); thref(n_points)];
    %         Ref = [Ref; ref];
    %     end
    % end

    % **************************** MPC ****************************

    solution = pred_control(x, A, B, C, Ref, Q, R, N, M, lb, ub);
    v = solution(1)+vref(i);
    w = solution(2)+wref(i);


    % for sim=1:2
    x(1) = x(1) + v*cos(x(3))*dt;
    x(2) = x(2) + v*sin(x(3))*dt;
    x(3) = x(3) + w*dt;
    % end
    vrobot(i) = v;
    wrobot(i) = w;

    % disp(A)
    % disp(B)
end

figure(1)
plot(xref, yref);
hold on;
plot(xrobot, yrobot);
legend('reference', 'real')
grid on;


figure(2)
plot(1:length(vrobot), vrobot);
hold on;
plot(1:length(wrobot), wrobot);
legend('v', 'w')
grid on;


figure(3)
plot(1:length(vref), vref);
hold on;
plot(1:length(wref), wref);
legend('vref', 'wref')
grid on;