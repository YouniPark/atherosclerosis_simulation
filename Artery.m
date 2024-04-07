classdef Artery
    % Model of artery (aortic or brachial)
    %   Detailed explanation goes here
    

    properties

        Initial_blood_volume {mustBeNumeric}
        Peripheral_resistance {mustBeNumeric}

        % aortic
        Length_aortic {mustBeNumeric}
        ELastic_modulus_aortic {mustBeNumeric}
        Wall_thickness_aortic {mustBeNumeric}
        Initial_diameter_aortic {mustBeNumeric}
        LDL_concentration_aortic {mustBeNumeric}
        

        % brachial
        Length_brachial {mustBeNumeric}
        ELastic_modulus_brachial {mustBeNumeric}
        Wall_thickness_brachial {mustBeNumeric}
        Initial_diameter_brachial {mustBeNumeric}
        LDL_concentration_brachial {mustBeNumeric}

        T_max {mustBeNumeric}
    end
    

    methods
        function obj = Artery(LDL_conc_aortic_, LDL_conc_brachial_)
           if nargin == 2
            
            obj.Initial_blood_volume = 500; 
            obj.Peripheral_resistance = 0.89*(0.22);
            obj.T_max = 0.8; 

            % Aortic:
            obj.LDL_concentration_aortic = LDL_conc_aortic_;

            obj.Length_aortic = 30; %mm
            obj.ELastic_modulus_aortic = 1.6; %MPa
            obj.Wall_thickness_aortic = 6; %mm
            obj.Initial_diameter_aortic = 36; %mm

            % Brachial:
            obj.LDL_concentration_brachial = LDL_conc_brachial_;

            obj.Length_brachial = 200; %mm
            obj.ELastic_modulus_brachial = 3.8; %MPa
            obj.Wall_thickness_brachial = 0.29; %mm
            obj.Initial_diameter_brachial = 4.3; %mm
           end
        end
        

        function [state_derivatives] = get_state_derivatives(obj, x, time)
            % calculate the time-varying blood flow
            blood_flow = obj.get_blood_flow(time);
            
            % calculate matrix A and B
            A = obj.matrix_A();
            B = obj.matrix_B();

            % implement the state equation x_dot = Ax + Bu
            state_derivatives = A * x + B * blood_flow;
            
            % Ensure state_derivatives is a column vector (correctly formatted)
            state_derivatives = state_derivatives(:);
                % Calculate the time-varying blood flo
        end
        

        function [A] = matrix_A(obj)
            % function to calculate matrix A

            % caluclate the resistances and compliances 
            % in aortic and brachial arteries
            resistances = get_resistance(obj);  % resistances = [Ra, Rb]
            compliances = get_compliances(obj); % compliances = [Ca, Cb]

            A = [1/(compliances(1) * resistances(2)), -1/(compliances(1) * resistances(2));
                -1/(compliances(2) * resistances(2)), 1/(compliances(2) * resistances(2)) - obj.Peripheral_resistance / resistances(2)];
           
        end


        function [B] = matrix_B(obj)
            % function to calculate matrix B
            
            % calculate the compliances
            compliances = get_compliances(obj); % compliances = [Ca, Cb]

            B = [1/compliances(1);
                0];

        end


        function [time_varying_blood_volume] = get_blood_volume(obj, time)
            % function to get blood volume wrt. time
            if time < 0.3
                time_varying_blood_volume = 500 * sin(pi * time / 0.3);
            else
                time_varying_blood_volume = 0;
            end
          
        end

        function [Blood_flow] = get_blood_flow(obj, time)
            %   Calculates finite-difference approximation of blood flow (blood volume derivative)
           
            % Input
            % time 
            
            % Output
            % finite-difference approximation of 
            % time derivative of time-varying blood volume (blood flow)
            
            if time < 0.3
                dt = 0.0001;
                forward_time = time + dt;
                backward_time = max(0, time - dt);
                forward = obj.get_blood_volume(forward_time);
                backward = obj.get_blood_volume(backward_time);
                Blood_flow = (forward - backward) / (2 * dt);
            else
                Blood_flow = 0;
            end
        end
        
        function [radius] = get_radius(obj)
            % function to get the radius of the artery wrt. time
            
            % Aortic:
            radius_a = obj.Initial_diameter_aortic / 2 * (1 - obj.LDL_concentration_aortic);

            % Brachial:
            radius_b = obj.Initial_diameter_brachial / 2 * (1 - obj.LDL_concentration_brachial);
            
            radius = [radius_a, radius_b];
           
        end

        function [compliances] = get_compliances(obj)
            % function to get compliances
            
            % calculate the current radius of the arteries
            radius = get_radius(obj);

            % Aortic:
            Ca = (3 * pi * (radius(1)) .^ 2) / (2 * obj.ELastic_modulus_aortic * obj.Wall_thickness_aortic);
            
            % Brachial:
            Cb = (3 * pi * (radius(2)) .^ 2) / (2 * obj.ELastic_modulus_brachial * obj.Wall_thickness_brachial);
            
            compliances = [Ca, Cb];
        end

        function [resistances] = get_resistance(obj)

            % function to get resistances
            
            % calculate the current radius of the arteries
            radius = get_radius(obj);

            % Aortic:
            resistance_a = (8 * obj.Length_aortic) / (pi * radius(1) .^ 4);

            % Brachial:
            resistance_b = (8 * obj.Length_brachial) / (pi * radius(2) .^ 4);

            resistances = [resistance_a, resistance_b];

        end


        function [time, y] = simulate(obj, total_time)
            % define the time span for simulation
            time_span = [0, total_time];

            % set the options to RelTol and AbsTol
            RelTol = 1e-6;
            AbsTol = 1e-8;
            options = odeset('RelTol', RelTol , 'AbsTol', AbsTol);

            % define initial blood pressure
            initial_blood_pressure =  [15, 15];

            % define the initial state
            [time, y] = ode45(@(t,x)obj.get_state_derivatives(x,t), time_span, initial_blood_pressure, options);

        end

        function [normalized_time] = get_normalized_time(obj, time)
            % DO LATER
            % Inputs
            % t: time
            
            % Output
            % time normalized to obj.Tmax (duration of each phase)
            
            normalized_time = mod(obj.Tmax, t);
        end
    end
end

