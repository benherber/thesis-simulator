classdef PathLoss
    %PATH_LOSS calcuations

    methods (Static)
        function res = amplitude_factor(distance, params)
            %AMPLITUDE_FACTOR due to friis free-space loss.

            arguments
                distance double
                params bherber.thesis.SimulationConstants
            end

            res = params.wavelen / (4 * pi * distance);
        end

        function res = power_factor(distance, params)
            %POWER_FACTOR due to friis free-space loss.

            arguments
                distance double
                params bherber.thesis.SimulationConstants
            end

            res = bherber.thesis.PathLoss.amplitude_factor(distance, params) .^ 2;
        end
    end
end