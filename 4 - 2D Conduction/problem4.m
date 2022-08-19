clear; clc; close all;

%% Parameters and Constants

k = 230;                               % Thermal conductivity [W/m K]
L.x = 0.10;                            % Length x [m]
L.y = 0.001 / 2;                       % Height y [m]
L.z = 1;                               % Width z [m]
P = 2 * L.z;                           % Fin perimeter [m]
Ac = L.z * 2*L.y;                      % Cross section area [m^2]
T_b = 100;                             % Base temperature [°C]
T_inf = 20;                            % External fluid temperature [°C]
eta = 0.6;                             % Fin efficiency
m = @(h) sqrt(h * P / (k * Ac));
h = fzero(@(h) eta - tanh(m(h) .* L.x) ./ (m(h) .* L.x), 100);
m = m(h);
M = sqrt(h * P * k * Ac) * (T_b - T_inf);
q_f = M * tanh(m * L.x)/2;             % Theoretical energy [W]

%% Discretization

Nx = 10;                               % Number of elements in x
Ny = Nx;                               % Number of elements y
N = Nx * Ny;                           % Total elements
grade = reshape(1:N, Nx, Ny)';
dx = L.x / (Nx - 0.5);                 % Distance between points in x   [m]
dy = L.y / (Ny - 1);                   % Distance between points in y   [m]
x = dx/2:dx:L.x;                       % Point coordinates in x [m]
y = 0:dy:L.y;                          % Point coordinates in y [m]
face.x = [0:dx:L.x-dx/2, L.x];         % Face coordinates in x  [m]
face.y = [0, dy/2:dy:L.y-dy/2, L.y];   % Face coordinates in y  [m]
delta.x = repmat(diff(face.x), 1, Ny); % Thickness of the volumes in x  [m]
delta.y = repelem(diff(face.y), Nx);   % Thickness of the volumes in y  [m]
area.W = delta.y .* L.z;               % Area of the west faces [m]
area.E = area.W;                       % Area of the east faces [m]
area.N = delta.x .* L.z;               % Area of the north faces [m]
area.S = area.N;                       % Area of the south faces   [m]
fronteira.W = grade(:, 1)';            % Indices for western border
fronteira.E = grade(:, Nx)';           % Indices for eastern border
fronteira.N = grade(Ny, :);            % Indices for northern border
fronteira.S = grade(1, :);             % Indices for southern border

% Boundary conditions
S_P = zeros(1, N);
S_P(fronteira.W) = S_P(fronteira.W) + -k .* area.W(fronteira.W) ./ (0.5 .* dx);
S_P(fronteira.E) = S_P(fronteira.E) + 0 .* area.E(fronteira.E);
S_P(fronteira.N) = S_P(fronteira.N) + -h .* area.N(fronteira.N);
S_P(fronteira.S) = S_P(fronteira.S) + 0 .* area.S(fronteira.S);
S_C = zeros(1, N);
S_C(fronteira.W) = S_C(fronteira.W) + k .* area.W(fronteira.W) .* T_b ./ (0.5 .* dx);
S_C(fronteira.E) = S_C(fronteira.E) + 0 .* area.E(fronteira.E);
S_C(fronteira.N) = S_C(fronteira.N) + h .* T_inf .* area.N(fronteira.N);
S_C(fronteira.S) = S_C(fronteira.S) + 0 .* area.S(fronteira.S);
b = S_C';

% Coefficients outside borders
for l = ["W", "E", "N", "S"]
    a.(l) = zeros(1, N);
    i = ~ismember(1:N, fronteira.(l));
    if l == "W" || l == "E"
        a.(l)(i) = k .* area.(l)(i) ./ dx;
    else
        a.(l)(i) = k .* area.(l)(i) ./ dy;
    end
    % Remove leading and trailing zeros
    a.(l) = a.(l)(find(a.(l), 1, 'first'):find(a.(l), 1, 'last'));
end

% Coefficients matrix (sparse)
A = diag(-a.E, 1) + diag(-a.W, -1) + diag(-a.N, Nx) + diag(-a.S, -Nx);
a.P = -sum(A) - S_P;
A = sparse(A + diag(a.P));

%% Solve for heat flow

T = A\b;
q_in = sum(k .* area.W(fronteira.W) .* (T_b - T(fronteira.W)') ./ (0.5 .* dx));
q_out = sum(h .* (T(fronteira.N)' - T_inf) .* area.N(fronteira.N));
balanco = q_in - q_out;
erro = abs(q_out - q_f) / q_f;

%% Visualization
T = reshape(T, [], Ny)';
Tm = mean(T);                   % Average temperatures

% T profile vs exact profile
figure('Name', "Temperature Profile")
X_an = linspace(0, L.x, 101);
T_an = T_inf + (T_b - T_inf) .* cosh(m .* (L.x - X_an)) ./ cosh(m .* L.x);
plot(x, Tm, 'ko', X_an, T_an, '-k', 'LineWidth', 1.25)
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12)
legend("Numeric", "Exact")
xlabel("Length (m)")
ylabel("Temperature (°C)")
box on; grid on;

% Isotherms
figure('Name', "Isotherms", 'Position', [25 50 1000 300])
[X, Y] = meshgrid([0 x], y .* 1e3);
Tc = [T_b.*ones(Ny, 1) T];
[C, hn] = contour(X, Y, Tc, 50:5:95, 'LineWidth', 1.5, 'LineColor', 'k');
clabel(C, hn, 'FontName', 'Times New Roman', 'FontSize', 12)
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12)
xlabel("Length (m)")
ylabel("Height (mm)")

% figure('Name', "Mesh")
% meshz(X, Y, Tc)
% xlabel("Length (m)", 'Rotation', -20)
% ylabel("Heigth (mm)", 'Rotation', 22)
% zlabel("Temperature (°C)")
% view(45, 30)

% Mesh visualization
% figure('Name', "Mesh", 'Position', [25 450 1000 300])
% [Xf,Yf] = ndgrid(face.x, face.y .* 1e3);
% plot(X, Y, 'ok', Xf, Yf, '-k', Xf', Yf', '-k', 'MarkerFaceColor', 'k');
% X = repmat(x, [1 Ny]);
% Y = repelem(y .* 1e3, Nx);
% str = {string(1:N)};
% text(X-0.001, Y+0.03, str{1}, 'FontName', 'Times New Roman')
% set(gca, 'FontName', 'Times New Roman', 'FontSize', 12)
% xlabel("Length (m)")
% ylabel("Heigth (mm)")
% ylim([0, L.y .* 1e3])
% box on;
