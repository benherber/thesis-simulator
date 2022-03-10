% --------------------------------------- %
% Frequency-Shift Keying Backscatter Tag  %
% Name   : FSKTag.m                       %
% Author : Benjamin Herber                %
% Date   : Spring 2021                    %
% --------------------------------------- %

classdef FSKTag < bherber.thesis.tags.Tag
    %FSK_TAG model
    %   Backscatter tag model of a passive Van Atta configuration using
    %   Frequency-Shift Keying to communicate

    %% Properties

    properties (GetAccess = public, SetAccess = private)
        f1 double {mustBePositive}
        f0 double {mustBePositive}
    end

    %% Public Instance Methods

    methods (Access = public)

        function this = FSKTag(x, y, z, mode, bits, sim_params, frequency_subband)
            %FSK_TAG constructor
            %   Create a tag a specific distance from a basestation in 3D
            %   space.
            
            arguments
                x double {mustBeScalarOrEmpty} % x-position of tag relative to basestation
                y double {mustBeScalarOrEmpty} % y-position of tag relative to basestation
                z double {mustBeScalarOrEmpty} % z-position of tag relative to basestation
                mode bherber.thesis.TagType % Modulation mode
                bits (1, :) {mustBeNonmissing, mustBeInRange(bits, 0, 1)} % bitstream
                sim_params bherber.thesis.SimulationConstants % Simulation parameters
                frequency_subband struct
            end

            this@bherber.thesis.tags.Tag(x, y, z, mode, bits, sim_params);
            this.f1 = frequency_subband.f1;
            this.f0 = frequency_subband.f0;

        end

% ----------------------------------------------------------------------- %
    
        function res = step(this)
            %STEP one simulation step using FSK modulation.
    
            if ~this.is_active
                    res = zeros(1, this.params.simstep_sz);
                    return
            end
    
            details = this.step_data();

            sq_one = square(2 * pi * this.f1 * details.time_slice);
            sq_zero = square(2 * pi * this.f0 * details.time_slice);
            
            all_ones = sq_one .* details.carrier_slice * ...
                this.vanatta_gain(this.params.num_elements, ...
                this.params.Fc, this.params.Fc + this.f1);
            all_zeros = sq_zero .* details.carrier_slice * ...
                this.vanatta_gain(this.params.num_elements, ...
                this.params.Fc, this.params.Fc + this.f0);

            res = zeros(size(details.data));
            for idx = 1:length(res)
                if details.data(idx) == 1
                    res(idx) = all_ones(idx);
                else
                    res(idx) = all_zeros(idx);
                end
            end

            res = res * (bherber.thesis.PathLoss.amplitude_factor(this.distance, this.params) .^ 2);
        end
    
% ----------------------------------------------------------------------- %
    
    function res_bits = demodulate(this, symb_num, signal)
            arguments
                this bherber.thesis.tags.FSKTag
                symb_num double {mustBePositive}
                signal (1, :)
            end

            time_slice = (0:(1 / this.params.Fs):(double(this.params.symb_sz) * (1 / this.params.Fs) - (1 / this.params.Fs))) ...
                + ((symb_num - 1) * double(this.params.symb_sz) * (1 / this.params.Fs));

            carrier_slice = this.params.amplitude * cos(2 * pi * this.params.Fc * time_slice);

            sq_one = square(2 * pi * this.f1 * time_slice);
            sq_zero = square(2 * pi * this.f0 * time_slice);
            
            all_ones = sq_one .* carrier_slice * ...
                this.vanatta_gain(this.params.num_elements, ...
                this.params.Fc, this.params.Fc + this.f1);
            all_zeros = sq_zero .* carrier_slice * ...
                this.vanatta_gain(this.params.num_elements, ...
                this.params.Fc, this.params.Fc + this.f0);

            % 1. Freq mix
            mixed_one = signal .* all_ones;
            mixed_zero = signal .* all_zeros;

            % 2. Filter
            filtered_one = lowpass(mixed_one, 2 * this.f1, this.params.Fs);
            filtered_zero = lowpass(mixed_zero, 2 * this.f0, this.params.Fs);

            % 3. Mean
            correlated_one = mean(filtered_one);
            correlated_zero = mean(filtered_zero);

            % 4. Combine Streams
            combined_streams = correlated_one - correlated_zero;

            % 5. Decide
            lambda = 0;
            if real(combined_streams) > lambda
                res_bits = 1;
            else
                res_bits = 0;
            end
        end
    end
end
