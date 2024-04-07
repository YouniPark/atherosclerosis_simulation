clear
clc

artery_model = Artery(0.9, 0);

% simulate for 0.3s
t = 0.3;
[time, state] = artery_model.simulate(t);

% calcultate the output y:

% calculate the aortic resistance
R = artery_model.get_resistance;
Ra = R(1);

% initialize a vector to hold the output
blood_pressure = zeros(length(time), 1);
blood_flow = zeros(length(time), 1);

% Loop over each time point
for i = 1:length(time)
    % Calculate blood flow at the current time point
    blood_flow(i) = artery_model.get_blood_flow(time(i));
    % Calculate 
    blood_pressure(i) = state(i,1) + Ra * blood_flow(i);
end

% Create data and 2-by-1 tiled chart layout
x = linspace(0,t,length(blood_flow));
tiledlayout(2,1)

% Plot blood flow
ax1 = nexttile;
plot(ax1,x,blood_flow)
title(ax1,'Blood Flow')
ylabel(ax1,'dV/dt')
xlabel(ax1, 'Time (seconds)')

% Plot blood_pressure
ax2 = nexttile;
plot(ax2,x,blood_pressure)
title(ax2,'Blood Pressure')
ylabel(ax2,'Pressure (mmHg)')
xlabel('Time (seconds)')
