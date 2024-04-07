clear
clc

artery_model = Artery(0.01, 0.01);

% simulate for 0.5s
t = 0.5;
[time, state] = artery_model.simulate(t);


% calcultate the output y:

% calculate the aortic resistance
R = artery_model.get_resistance;
Ra = R(1);

% initialize a vector to hold the output
y = zeros(length(time), 1);

% Loop over each time point
for i = 1:length(time)
    % Calculate blood flow at the current time point
    Qt = artery_model.get_blood_flow(time(i));
    % Calculate y
    y(i) = state(i,1) + Ra * Qt;
end



% Plot
figure()
LineWidth = 1.5;
FontSize = 12;

plot(time, y, 'r', 'LineWidth', LineWidth)
legend('blood pressure') % note the order here
xlabel('Time (seconds)')
ylabel('Pressure (mmHg)')
set(gca, 'FontSize', FontSize)