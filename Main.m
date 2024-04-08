clear
clc

% new radius = initial radius(1-factor)
artery_model = Artery(0, 0.8); 

% simulate for 0.8s
t = 0.8;
% normalized_time = artery_model.get_normalized_time(t);
[time, state] = artery_model.simulate(t);

blood_volume = artery_model.get_blood_volume(0.29);

% calculate the aortic resistance
R = artery_model.get_resistance;
Ra = R(1);
Rb = R(2);

% get compliance for aortic
compliances = artery_model.get_compliances;
Ca = compliances(1);
Cb = compliances(2);

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

%% Unit tests:
clear
clc

results = runtests('Arterytests');
    disp(results)

%% Sensitivity Checks:

% test on sensitivity to aortic percent reduction
percentReductions = linspace(0, 0.5, 10); % 0% to 50% reduction
bloodPressureResponses = zeros(length(percentReductions), 1);

for i = 1:length(percentReductions)
    model = Artery(percentReductions(i), 0.8); % Keeping brachial reduction constant
    [~, state] = model.simulate(0.8); 
    bloodPressureResponses(i) = mean(state(:,1)); % Average aortic pressure
end

figure;
plot(percentReductions, bloodPressureResponses);
title('Sensitivity of Blood Pressure to Aortic Area Reduction');
xlabel('Aortic Area Reduction (%)');
ylabel('Average Aortic Blood Pressure (mmHg)');


% test on sensitivity to elastic moduli
elasticModuli = linspace(1, 4, 5); % range from 1 MPa to 4 MPa
averageBloodPressureAortic = zeros(length(elasticModuli), 1);

for i = 1:length(elasticModuli)
    model = Artery(0.2, 0.8); 
    model.ELastic_modulus_aortic = elasticModuli(i);
    model.ELastic_modulus_brachial = elasticModuli(i); % Assuming the same change
    [~, state] = model.simulate(0.8);
    averageBloodPressureAortic(i) = mean(state(:,1)); % Average aortic pressure
end

figure;
plot(elasticModuli, averageBloodPressureAortic);
title('Sensitivity of Blood Pressure to Elastic Modulus');
xlabel('Elastic Modulus (MPa)');
ylabel('Average Aortic Blood Pressure (mmHg)');


% test on sensitivity to wall thickness
wallThicknesses = linspace(0.1, 10, 5); % range from 0.1 mm to 10 mm
averageBloodPressureAortic = zeros(length(wallThicknesses), 1);

for i = 1:length(wallThicknesses)
    model = Artery(0.2, 0.8); 
    model.Wall_thickness_aortic = wallThicknesses(i);
    model.Wall_thickness_brachial = wallThicknesses(i); % Assuming the same change for simplicity
    [~, state] = model.simulate(0.8);
    averageBloodPressureAortic(i) = mean(state(:,1)); % Average aortic pressure
end

figure;
plot(wallThicknesses, averageBloodPressureBrachial);
title('Sensitivity of Blood Pressure to Wall Thickness');
xlabel('Wall Thickness (mm)');
ylabel('Average Aortic Blood Pressure (mmHg)');


% test on sensitivity to peripheral resistances
peripheralResistances = linspace(0.1, 1.5, 5); 
averageSystemicPressure = zeros(length(peripheralResistances), 1);

for i = 1:length(peripheralResistances)
    model = Artery(0.2, 0.8); 
    model.Peripheral_resistance = peripheralResistances(i);
    [~, state] = model.simulate(0.8);
    % Considering systemic pressure as an average of aortic and brachial pressures
    averageSystemicPressure(i) = mean(state(:,1) + state(:,2)) / 2;
end

figure;
plot(peripheralResistances, averageSystemicPressure);
title('Sensitivity of Systemic Blood Pressure to Peripheral Resistance');
xlabel('Peripheral Resistance');
ylabel('Average Systemic Blood Pressure (mmHg)');