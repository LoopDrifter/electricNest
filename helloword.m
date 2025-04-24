%code1
bmax = 1;           
freq = 50;         
w = 2*pi*freq;       


t = 0:1/5000:1/50;
Baa = sin(w*t) .* (cos(0) + j*sin(0));
Bbb = sin(w*t-2*pi/3) .* (cos(2*pi/3) + j*sin(2*pi/3));
Bcc = sin(w*t+2*pi/3) .* (cos(-2*pi/3) + j*sin(-2*pi/3));

Bnet = Baa + Bbb + Bcc;

circle = 1.5 * (cos(w*t) + j*sin(w*t));

for ii = 1:length(t)
    plot(circle, 'g');
    hold on;

    plot([0 real(Baa(ii))], [0 imag(Baa(ii))], 'k', 'LineWidth', 2);
    plot([0 real(Bbb(ii))], [0 imag(Bbb(ii))], 'b', 'LineWidth', 2);
    plot([0 real(Bcc(ii))], [0 imag(Bcc(ii))], 'm', 'LineWidth', 2);
    plot([0 real(Bnet(ii))], [0 imag(Bnet(ii))], 'r', 'LineWidth', 3);
    
    axis square;
    axis([-2 2 -2 2]);
    drawnow;
    hold off;
end



%code 2



i_a = (0:1:20) * 3;


e_a = 277.0;  
x_s = 1.0;    

pf_values = [0.2, 0.4, 0.6, 0.8];
theta_lagging = acos(pf_values);   
theta_leading = acos(-pf_values);  

colors = ['r', 'g', 'b', 'm'];

figure;
hold on;


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

xlabel('Line Current (A)', 'Fontweight', 'Bold');
ylabel('Terminal Voltage (V)', 'Fontweight', 'Bold');
title('Terminal Characteristics for Various Power Factors', 'Fontweight', 'Bold');
grid on;
axis([0 60 400 550]);
legend show;
hold off;




%cide 3
clc;
clear;


V = 460;              
f = 60;               
p = 4;                 
Ns = 120 * f / p;      
ws = 2 * pi * Ns / 60; 


XR = 0.8;            


RR_base = 0.4;       
RR_half = RR_base / 2;
RR_double = RR_base * 2;


s = linspace(0.001, 1, 1000);


T_base = (3 * V^2 .* (RR_base ./ s)) ./ (ws * ((RR_base ./ s).^2 + XR^2));
T_half = (3 * V^2 .* (RR_half ./ s)) ./ (ws * ((RR_half ./ s).^2 + XR^2));
T_double = (3 * V^2 .* (RR_double ./ s)) ./ (ws * ((RR_double ./ s).^2 + XR^2));

N = (1 - s) * Ns;


figure;
plot(N, T_base, 'b', 'LineWidth', 2); hold on;
plot(N, T_half, 'r--', 'LineWidth', 2);
plot(N, T_double, 'g-.', 'LineWidth', 2);

xlabel('Rotor Speed (RPM)');
ylabel('Torque (Nm)');
title('Torque-Speed Characteristics of an Induction Motor');
legend('R_R', '0.5 R_R', '2 R_R', 'Location', 'Best');
grid on;
