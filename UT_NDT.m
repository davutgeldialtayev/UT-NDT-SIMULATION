clc;
clear;
close all;

%% Material Properties
density = 7800;        % Density of steel in kg/m^3
speed_of_sound = 5000; % Speed of ultrasonic wave in m/s
poisson_ratio = 0.3;   % Poisson's ratio for steel
young_modulus = 210e9; % Young's modulus for steel (in Pascals)
ultrasonic_freq = 1e6; % Frequency of ultrasonic waves in Hz

%% Simulation Parameters
material_thickness = 0.01;   % Material thickness in meters (10 mm)
flaw_depth = 0.005;          % Depth of internal flaw in meters (5 mm)
flaw_size = 0.002;           % Size of the flaw (diameter) in meters (2 mm)
distance_to_flaw = flaw_depth; % Distance to flaw

% Initial Echo Time (without flaw)
echo_time_no_flaw = 2 * material_thickness / speed_of_sound; % Time for wave to return from surface
echo_time_with_flaw = 2 * (material_thickness + flaw_depth) / speed_of_sound; % Time for wave to return from flaw

% Number of steps for thickness variation
num_steps = 10;
thickness_values = linspace(0.005, 0.02, num_steps); % Thickness from 5 mm to 20 mm

%% Preallocate Echo Time Array
echo_time_no_flaw_arr = zeros(num_steps, 1);
echo_time_with_flaw_arr = zeros(num_steps, 1);

% Simulate Echo Time for Different Material Thicknesses
for i = 1:num_steps
    current_thickness = thickness_values(i);
    echo_time_no_flaw_arr(i) = 2 * current_thickness / speed_of_sound;  % Echo time for no flaw
    echo_time_with_flaw_arr(i) = 2 * (current_thickness + flaw_depth) / speed_of_sound; % Echo time for flaw
end

%% 3D Visualization of Ultrasonic Wave Propagation through Material

% Create a grid for 3D plot
[X, Y] = meshgrid(linspace(-material_thickness, material_thickness, 100), ...
                  linspace(0, material_thickness + flaw_depth, 100));

% Simulate wave propagation (this is a simple representation for visualization)
Z_no_flaw = exp(-((X.^2 + Y.^2) / (2 * (material_thickness / 4)^2))); % Wave without flaw
Z_with_flaw = exp(-((X.^2 + (Y - flaw_depth).^2) / (2 * (material_thickness / 4)^2))); % Wave with flaw

% 3D Surface plot for wave propagation (no flaw)
figure;
surf(X, Y, Z_no_flaw, 'EdgeColor', 'none');
colormap jet;
title('3D Ultrasonic Wave Propagation Without Flaw');
xlabel('Material Thickness (m)');
ylabel('Depth (m)');
zlabel('Wave Amplitude');
colorbar;
view(30, 30);
grid on;

% 3D Surface plot for wave propagation (with flaw)
figure;
surf(X, Y, Z_with_flaw, 'EdgeColor', 'none');
colormap jet;
title('3D Ultrasonic Wave Propagation With Flaw');
xlabel('Material Thickness (m)');
ylabel('Depth (m)');
zlabel('Wave Amplitude');
colorbar;
view(30, 30);
grid on;

%% Echo Time vs. Material Thickness (2D Plot)
figure;
plot(thickness_values * 1000, echo_time_no_flaw_arr * 1000, 'b-', 'LineWidth', 2); % in milliseconds
hold on;
plot(thickness_values * 1000, echo_time_with_flaw_arr * 1000, 'r--', 'LineWidth', 2); % with flaw
xlabel('Material Thickness (mm)');
ylabel('Echo Time (ms)');
legend('No Flaw', 'With Flaw');
title('Echo Time vs. Material Thickness');
grid on;

%% Display Results for a 10mm Thickness
fprintf('For a material thickness of 10 mm:\n');
fprintf('Echo time without flaw: %.2f ms\n', echo_time_no_flaw_arr(5) * 1000);
fprintf('Echo time with flaw: %.2f ms\n', echo_time_with_flaw_arr(5) * 1000);

%% 3D Waveform Visualization for No Flaw and With Flaw
time_vector = linspace(0, echo_time_with_flaw_arr(5) * 1.2, 500); % Reduce points to manage memory
no_flaw_waveform = sin(2 * pi * ultrasonic_freq * time_vector); % No flaw waveform
with_flaw_waveform = sin(2 * pi * ultrasonic_freq * time_vector); % Flaw waveform

% Adjust the time vector to match the echo times and visualize the waveform delay
delay_time = echo_time_with_flaw_arr(5) * 1.2 / 2; % Half of the echo time for reasonable delay

% Apply delay to waveforms
no_flaw_waveform_delay = circshift(no_flaw_waveform, round(delay_time * length(time_vector) / (time_vector(end) - time_vector(1))));
with_flaw_waveform_delay = circshift(with_flaw_waveform, round(delay_time * length(time_vector) / (time_vector(end) - time_vector(1))));

% Create a 3D plot for the waveforms
figure;
subplot(2, 1, 1);
plot3(time_vector, no_flaw_waveform_delay, zeros(size(time_vector)), 'b-', 'LineWidth', 2);
title('3D Ultrasonic Waveform Without Flaw');
xlabel('Time (s)');
ylabel('Amplitude');
zlabel('Waveform ID');
grid on;

subplot(2, 1, 2);
plot3(time_vector, with_flaw_waveform_delay, ones(size(time_vector)), 'r--', 'LineWidth', 2);
title('3D Ultrasonic Waveform With Flaw');
xlabel('Time (s)');
ylabel('Amplitude');
zlabel('Waveform ID');
grid on;

%% 3D Flaw Detection Visualization
% Visualize the flaw inside the material in 3D space
[X_flaw, Y_flaw] = meshgrid(linspace(0, material_thickness, 50), linspace(0, flaw_depth, 50));
Z_flaw = zeros(size(X_flaw)); % Represent the flaw as a flat region in material

figure;
mesh(X_flaw, Y_flaw, Z_flaw, 'EdgeColor', 'r', 'LineWidth', 3);
hold on;
mesh(X, Y, Z_no_flaw, 'EdgeColor', 'b');
xlabel('Material Thickness (m)');
ylabel('Flaw Depth (m)');
zlabel('Flaw Location');
title('3D Visualization of Flaw Inside Material');
legend('Flaw', 'Material');
grid on;



