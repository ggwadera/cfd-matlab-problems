clear;

%% Parameters and Constants

k  = 1;                          % Conductivity - W/(m*K)
cp = 1000;                       % Specific Heat - J/(kg*K)
rho = 3000;                      % Specific Mass - kg/m^3
L  = 0.2;                        % Lenght - m
Bi = 1.3;                        % Biot Number
h  = Bi * k / L;                 % Convective Coefficient - W/(m^2*K)
T_inf = 200;                     % Fluid temperate - °C
Ti = 25;                         % Initial temperature - °C
tf = 200 * 3600;                 % Final time - s
N  = 100;                        % Number of volumes
dx = L / (N - 1);                % Volume length - m
dt = 0.45 * dx^2 * rho * cp / k; % Time step - s
x  = 0:dx:L;                     % Coordinates vector
t  = 0:dt:tf;                    % Time vector

%% Solving

T = ones(length(t), N) .* Ti;    % Temperature matix
a = k / dx;                      % Coefficients aw and ae
ap = rho * cp * dx / dt;         % Coefficients ap and ap_0

% Temperature profile calculation by time stepping
for j = 1:length(t)-1
    % Temperature at left surface
    T(j+1, 1) = (a * T(j, 2) + (ap/2 - a) * T(j, 1)) / (ap/2);
    
    % Temperatures on the interior
    for i = 2:N-1
        T(j+1, i) = (a * (T(j, i-1) + T(j, i+1)) + ...
                    (ap - 2 * a) * T(j,i)) / ap;
    end
    
    % Temperature at right surface
    T(j+1, N) = (a * T(j, N-1) + (ap/2 - h - a) * ...
                T(j, N) + h * T_inf) / (ap/2);
end

%% Results

line = ["-*", "-o", "-d", "-s", "-+"]; l = 0;
figure(); hold on; grid on;
for i = 1:floor(length(t)/4):length(t)
    l = l + 1;
    plot(x, T(i, :), 'k' + line(l), 'LineWidth', 1.5, ...
        'MarkerIndices', 1:N/10:N, ...
        'DisplayName', "t = " + round(t(i)/3600, 2) + " h")
end
xlabel("Length (m)");
ylabel("Temperature (°C)");
ylim([20, 200])
% yticks(20:10:200)
legend('Location', 'southeast');