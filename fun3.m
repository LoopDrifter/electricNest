clc;
clear;

% Given parameters
V = 460;               % Supply voltage in Volts
f = 60;                % Frequency in Hz
p = 4;                 % Number of poles
Ns = 120 * f / p;      % Synchronous speed in RPM
ws = 2 * pi * Ns / 60; % Synchronous speed in rad/s

% Rotor reactance (constant)
XR = 0.8;              % Ohms

% Rotor resistance values
RR_base = 0.4;         % Chosen value
RR_half = RR_base / 2;
RR_double = RR_base * 2;

% Slip range
s = linspace(0.001, 1, 1000);

% Torque for each resistance
T_base = (3 * V^2 .* (RR_base ./ s)) ./ (ws * ((RR_base ./ s).^2 + XR^2));
T_half = (3 * V^2 .* (RR_half ./ s)) ./ (ws * ((RR_half ./ s).^2 + XR^2));
T_double = (3 * V^2 .* (RR_double ./ s)) ./ (ws * ((RR_double ./ s).^2 + XR^2));

% Speed vector
N = (1 - s) * Ns;

% Plotting
figure;
plot(N, T_base, 'b', 'LineWidth', 2); hold on;
plot(N, T_half, 'r--', 'LineWidth', 2);
plot(N, T_double, 'g-.', 'LineWidth', 2);

xlabel('Rotor Speed (RPM)');
ylabel('Torque (Nm)');
title('Torque-Speed Characteristics of an Induction Motor');
legend('R_R', '0.5 R_R', '2 R_R', 'Location', 'Best');
grid on;