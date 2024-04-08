classdef Arterytests < matlab.unittest.TestCase
    % Test class for Artery.m
    
    properties
        Artery
    end

    methods(TestMethodSetup)
        function createArtery(testCase)
            testCase.Artery = Artery(0, 0.8); % Initialize with randomly selected percent initial radius
        end
    end

    methods(TestMethodTeardown)
    end

    methods(Test)
        function testBloodVolumeDecreasePhase(testCase)
            % Verify blood volume decreases after a certain time
           
            time = 0.4; % After t = 0.3s
            vol = testCase.Artery.get_blood_volume(time);
            testCase.assertLessThan(vol, testCase.Artery.Initial_blood_volume, 'Blood volume should decrease in the latter phase');
        end

        function testBloodFlowSignChange(testCase)
            % Test blood flow changes sign (indicating direction change)

            time1 = 0.15;
            time2 = 0.40;
            flow1 = get_blood_flow(testCase.Artery, 0.1);
            flow2 = get_blood_flow(testCase.Artery, 0.3);
            display(flow1);
            display(flow2);
            % checking for sign change
            hasSignChanged = sign(flow1) ~= sign(flow2);
            testCase.assertTrue(hasSignChanged, 'Blood flow should change direction');
        end
        
        function testStateDerivativesSize(testCase)
            % Ensure get_state_derivatives returns a vector of correct size

            x = [80; 85]; % Initial state
            time = 0.1;
            state_derivatives = testCase.Artery.get_state_derivatives(x, time);
            testCase.assertSize(state_derivatives, [2, 1], 'State derivatives vector should be of size 2x1');
        end

        function testCompliancesCalculation(testCase)
            % Test compliances calculation
            model_compliances = testCase.Artery.get_compliances();
            expected_compliances = [636.17,3.16];  % obtained by hand calculation
            testCase.verifyEqual(model_compliances, expected_compliances, 'AbsTol', 0.01);
            
        end

        function testResistancesCalculation(testCase)
            % Test resistances calculation
            model_resistances = testCase.Artery.get_resistance();
            expected_resistances = [0,1082.62];  % obtained by hand calculation
            testCase.verifyEqual(model_resistances, expected_resistances, 'AbsTol', 0.01);
        end

    end
    
end