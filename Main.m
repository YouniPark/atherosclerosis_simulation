clear
clc

% new radius = initial radius(1-factor)
artery_model = Artery(0, 0.8); 

% simulate for 0.3s
t = 0.8;
% normalized_time = artery_model.get_normalized_time(t);
[time, state] = artery_model.simulate(t);

blood_volume = artery_model.get_blood_volume(0.29);

% calculate the aortic resistance
R = artery_model.get_resistance;
Ra = R(1)
Rb = R(2)

% get compliance for aortic
compliances = artery_model.get_compliances;
Ca = compliances(1)
Cb = compliances(2)

% initialize a vector to hold the output

blood_pressure = zeros(length(time), 1);
blood_pressure_a = zeros(length(time), 1);
blood_pressure_b = zeros(length(time), 1);
blood_flow = zeros(length(time), 1);

% Loop over each time point
for i = 1:length(time)
    % Calculate blood flow at the current time point
    blood_flow(i) = artery_model.get_blood_flow(time(i));
    % Calculate 
    % blood_pressure(i) = ( (Ra + Ca) * blood_flow(i) + Ca / Rb * (state(i, 1) - state(i, 2)) );
    blood_pressure(i) = state(i, 1) + 1/3 * state(i, 2) + Ra * blood_flow(i);
    blood_pressure_a(i) = state(i,1);
    blood_pressure_b(i) = state(i,2);
end

% Create data and 2-by-1 tiled chart layout
x = linspace(0,t,length(blood_flow));
tiledlayout(4,1)

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

% Plot blood_pressure in aortic arch
ax3 = nexttile;
plot(ax3,x,blood_pressure_a)
title(ax3,'Blood Pressure in Aortic Arch')
ylabel(ax3,'Pressure (mmHg)')
xlabel('Time (seconds)')

% Plot blood_pressure in brachial artery
ax4 = nexttile;
plot(ax4,x,blood_pressure_b)
title(ax4,'Blood Pressure in Brachial Artery')
ylabel(ax4,'Pressure (mmHg)')
xlabel('Time (seconds)')
