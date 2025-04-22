% M-file: term_char_all.m
% M-file to plot the terminal characteristics of the generator
% for various power factors (both leading and lagging)

% Initialize the current amplitudes (21 values in the range 0-60 A)
i_a = (0:1:20) * 3;

% Initialize other parameters
e_a = 277.0;  % Internal generated voltage
x_s = 1.0;    % Synchronous reactance

% Define power factor angles (in radians) for lagging and leading cases
pf_values = [0.2, 0.4, 0.6, 0.8];
theta_lagging = acos(pf_values);   % Lagging power factor (current lags voltage)
theta_leading = acos(-pf_values);  % Corrected Leading power factor calculation

% Colors for different plots
colors = ['r', 'g', 'b', 'm'];

figure;
hold on;

% Plot terminal characteristics for lagging power factors
for idx = 1:length(pf_values)
    v_phase = zeros(1, 21);
    theta = theta_lagging(idx);
    
    for ii = 1:21
        v_phase(ii) = sqrt(e_a^2 - (x_s * i_a(ii) * sin(theta))^2) ...
                    - (x_s * i_a(ii) * cos(theta));
    end
    
    v_t = v_phase * sqrt(3);
    plot(i_a, v_t, 'Color', colors(idx), 'Linewidth', 2.0, 'DisplayName', ...
        sprintf('Lagging PF = %.1f', pf_values(idx)));
end

% Plot terminal characteristics for leading power factors
for idx = 1:length(pf_values)
    v_phase = zeros(1, 21);
    theta = theta_leading(idx);
    
    for ii = 1:21
        v_phase(ii) = sqrt(e_a^2 - (x_s * i_a(ii) * sin(theta))^2) ...
                    - (x_s * i_a(ii) * cos(theta));
    end
    
    v_t = v_phase * sqrt(3);
    plot(i_a, v_t, '--', 'Color', colors(idx), 'Linewidth', 2.0, 'DisplayName', ...
        sprintf('Leading PF = %.1f', pf_values(idx)));
end

% Add labels, title, and legend
xlabel('Line Current (A)', 'Fontweight', 'Bold');
ylabel('Terminal Voltage (V)', 'Fontweight', 'Bold');
title('Terminal Characteristics for Various Power Factors', 'Fontweight', 'Bold');
grid on;
axis([0 60 400 550]);
legend show;
hold off;