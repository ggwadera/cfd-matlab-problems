%% Parameters and Constants

ki = [10, 1];                   % Thermal Conductivities - W/(m*K)
Li = [0.069, 0.031];            % Lenghts - m
Ni = [3, 3];                    % Number of volumes
h   = 100;                      % Convective coefficient - W/(m^2 * K)
T_inf = 40;                     % Fluid temperature - °C
q_presc = 6000;                 % Prescribed heat flux - W/m^2
Area = 1;                       % Surface area - m^2
exact = [327.4, 100];           % Exact solution values

%% Data initialization

N  = sum(Ni);                        % Total volume numbers
dx = repelem(Li ./ (Ni - 0.5), Ni);  % Width of each volume
k  = repelem(ki, Ni);                % k of each point
b1 = q_presc;                        % boundary condition at the start
bn = h * Area * T_inf;               % boundary condition at the end
b = [b1; zeros(N-2, 1); bn];         % boundary conditions array

%% Discretization

% Interface conductivity by thermal resistance method
fe = dx(2:N) ./ (dx(1:N-1) + dx(2:N));
fw = dx(1:N-1) ./ (dx(1:N-1) + dx(2:N));
ke = 1 ./ ((1 - fe) ./ k(1:N-1) + fe ./ k(2:N));
kw = 1 ./ ((1 - fw) ./ k(2:N) + fw ./ k(1:N-1));
dx_PE = 0.5 .* (dx(1:N-1) + dx(2:N));
dx_PW = 0.5 .* (dx(1:N-1) + dx(2:N));
ae = ke .* Area ./ dx_PE;
aw = kw .* Area ./ dx_PW;
ap = [ae, h*Area] + [0, aw];
A = diag(-ae, 1) + diag(ap) + diag(-aw, -1);

% Interface conductivity by linear interpolation
ke_lin = fe .* k(1:N-1) + (1 - fe) .* k(2:N);
kw_lin = fw .* k(2:N) + (1 - fw) .* k(1:N-1);
ae_lin = ke_lin .* Area ./ dx_PE;
aw_lin = kw_lin .* Area ./ dx_PW;
ap_lin = [ae_lin, h*Area] + [0, aw_lin];
A_lin = diag(-ae_lin, 1) + diag(ap_lin) + diag(-aw_lin, -1);

% Solve linear system of equations
T = A \ b;
T_lin = A_lin \ b;
erro = [T(1), T(N); T_lin(1), T_lin(N)] - exact;

%% Plot temperature profile

% Coordinates in mm
x = cumsum([0, dx_PE]) .* 1000;

figure(); hold on; grid on;
plot(x, T, 'ko', 'MarkerSize', 7, 'MarkerFaceColor', 'black')
plot(x, T_lin, 'k*', 'MarkerSize', 12, 'MarkerFaceColor', 'black')
plot([0, 69, 100], [327.4, 286, 100], '-k', 'LineWidth', 1.5)
xlabel('Length (mm)')
ylabel('Temperature (°C)')
legend('Resistances', 'Linear', 'Exact', 'Location', 'southwest');
