clear;
% clc
%% Parameters and Constants

k   = 230;                        % Conductivity - W/(m*K)
rho = 2700;
cp  = 900;
R0  = 50e-3;                      % Length - m
R   = 55e-3;
L   = R - R0;
N   = 100;                        % Number of volumes in A and B
h   = 20000;                      % Convective coefficient - W/(m^2 * K)
T_inf = 80;                       % Fluid temperate - °C
T_sup = @(t) 150 + 50 .* cos(50 .* pi .* t);
T0    = 80;
H     = 1;

%% Discretization and initialization

dx = L / (N - 1);     % Distance between points
r  = R0:dx:R;         % Point coordinates
dt = 0.02 / 100;
tf = dt * 1000;
T(1:N, 1) = T0;

%% Solving

A = zeros(N, N);
A(1, 1) = 1;
b = zeros(1, N);
b_N = h * T_inf * 2 * pi * R * H;

for t = 0:dt:tf
    b(1) = T_sup(t);
    for i = 2:N
        if i == N
            Area_w = 2 * pi * (R - dx/2) * H;
            Vol = pi * (R^2 - (R-dx/2)^2) * H;
            Mp = rho * Vol;
            aw = k * Area_w / dx;
            ae = h * 2 * pi * R * H;
            ap0 = Mp * cp / dt;
            ap = ae + aw + ap0;
            A(i, [i-1, i]) = [-aw, ap];
            b(i) = b_N + ap0 * T(i);
        else
            Area_w = 2 * pi * (r(i) - dx/2) * H;
            Area_e = 2 * pi * (r(i) + dx/2) * H;
            Vol = pi * ((r(i)+dx/2)^2 - (r(i)-dx/2)^2) * H;
            Mp = rho * Vol;
            aw = k * Area_w / dx;
            ae = k * Area_e / dx;
            ap0 = Mp * cp / dt;
            ap = ae + aw + ap0;
            A(i, [i-1, i, i+1]) = [-aw, ap, -ae];
            b(i) = ap0 * T(i);
        end
    end

    % Solve linear equations
    T = A \ b';
%     plot(t, T(5), 'k.')
%     plot(t, T(N), 'b.')
end

%% Plot temperate profile

% % Point coordinates in mm
x = cumsum([0, dx_PE]) .* 1000;

figure(); hold on; grid on;
plot(x, T, 'ko', 'MarkerSize', 7, 'MarkerFaceColor', 'black')
plot(x, T_lin, 'k*', 'MarkerSize', 12, 'MarkerFaceColor', 'black')
plot([0, 69, 100], [327.4, 286, 100], '-k', 'LineWidth', 1.5)
xlabel('Length (mm)')
ylabel('Temperature (°C)')
legend('Resistances', 'Linear', 'Exact', 'Location', 'southwest');