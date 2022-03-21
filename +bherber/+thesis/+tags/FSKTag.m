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

        function this = FSKTag(x, y, z, mode, bits, sim_params, frequency_subband, options)
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
                options.snr_db = NaN;
                options.complex_noise = false;
            end

            this@bherber.thesis.tags.Tag(x, y, z, mode, bits, sim_params, ...
                snr_db=options.snr_db, complex_noise=options.complex_noise);
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

            pathloss = bherber.thesis.PathLoss.amplitude_factor(this.distance, this.params);
            gain = this.vanatta_gain(this.params.num_elements, this.params.Fc, this.params.Fc);

            sq_one = square(2 * pi * this.f1 * details.time_slice);
            sq_zero = square(2 * pi * this.f0 * details.time_slice);
            
            all_ones = (sq_one .* details.carrier_slice) .* gain;
            all_zeros = (sq_zero .* details.carrier_slice) .* gain;

            res = zeros(size(details.data));
            for idx = 1:length(res)
                if details.data(idx) == 1
                    res(idx) = all_ones(idx);
                else
                    res(idx) = all_zeros(idx);
                end
            end

            if ~isnan(this.snr_db)
                linear_snr = 10 ^ (this.snr_db / 10);
                expected_amplitude = this.params.amplitude * pathloss * gain;
                linear_carrier_power = (expected_amplitude ^ 2) / 2; % ASSUMPTION: Uniform Random Data
                signal_noise = bherber.thesis.channel_models.AWGNChannelModel.noise(...
                    length(details.time_slice), linear_carrier_power, linear_snr, this.complex_noise);
                res = res + signal_noise;
            end

            res = res * pathloss;
        end

% ----------------------------------------------------------------------- %
    
        function point = constellation_point(this, symb_num, signal)
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
            
            all_ones = sq_one .* carrier_slice;
            all_zeros = sq_zero .* carrier_slice;

            % 1. Freq mix
            mixed_one = signal .* all_ones;
            mixed_zero = signal .* all_zeros;

            % 2. Filter
            filtered_one = lowpass(mixed_one, 2 * this.f1, this.params.Fs);
            filtered_zero = lowpass(mixed_zero, 2 * this.f0, this.params.Fs);

            % 3. Mean
            correlated_one = trapz(filtered_one);
            correlated_zero = trapz(filtered_zero);

            % 4. Combine Streams
            combined_streams = abs(correlated_one) - abs(correlated_zero);
            point = combined_streams;
        end
    
% ----------------------------------------------------------------------- %
    
    function res_bits = demodulate(this, symb_num, signal)
            arguments
                this bherber.thesis.tags.FSKTag
                symb_num double {mustBePositive}
                signal (1, :)
            end

            combined_streams = this.constellation_point(symb_num, signal);

%             gca;
%             hold on
%             ramp = [];
%             for idx = 1:100:(length(time_slice) - 1)
%                 ramp = [ramp, (1 / this.params.Fs) * ...
%                     (abs(trapz(filtered_one(1:idx))) - abs(trapz(filtered_zero(1:idx))))];
%             end
%             if isnan(this.snr_db)
%                 plot(gca, ramp, "-b");
%             else
%                 plot(gca, ramp, "-r");
%             end

            % 5. Decide
            lambda = 0;
            if combined_streams > lambda
                res_bits = 1;
            else
                res_bits = 0;
            end
        end
    end
end
