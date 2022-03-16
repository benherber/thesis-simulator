% ---------------------------------- %
% Frequency-Hopping Backscatter Tag  %
% Name   : FreqHopTag.m              %
% Author : Benjamin Herber           %
% Date   : Spring 2021               %
% ---------------------------------- %

classdef FreqHopTag < bherber.thesis.tags.Tag
    %FREQ_HOP_TAG model
    %   Backscatter tag model of a passive Van Atta configuration using
    %   Frequency-Hopping to communicate.

    %% Properties

    properties (GetAccess = public, SetAccess = private)
        hop_pattern (1, :) double {mustBePositive}
        freq_subbands (1, :) struct
    end

    %% Public Instance Methods

    methods (Access = public)

        function this = FreqHopTag(x, y, z, mode, bits, sim_params, frequency_subbands, hop_pattern)
            %FREQ_HOP_TAG constructor
            %   Create a tag a specific distance from a basestation in 3D
            %   space.
            
            arguments
                x double {mustBeScalarOrEmpty} % x-position of tag relative to basestation
                y double {mustBeScalarOrEmpty} % y-position of tag relative to basestation
                z double {mustBeScalarOrEmpty} % z-position of tag relative to basestation
                mode bherber.thesis.TagType % Modulation mode
                bits (1, :) {mustBeNonmissing, mustBeInRange(bits, 0, 1)} % bitstream
                sim_params bherber.thesis.SimulationConstants % Simulation parameters
                frequency_subbands (1, :) struct
                hop_pattern (1, :)
            end

            this@bherber.thesis.tags.Tag(x, y, z, mode, bits, sim_params);
            this.freq_subbands = frequency_subbands;
            this.hop_pattern = hop_pattern;

        end

% ----------------------------------------------------------------------- %
    
        function res = step(this)
            %STEP one simulation step using FSK modulation.
    
            if ~this.is_active
                res = zeros(1, this.params.simstep_sz);
                return
            end
    
            details = this.step_data();
            delay_samples = details.delay_samples;
            pattern_len = length(this.hop_pattern);

            pattern_idx = mod(details.curr_symb, pattern_len);
            if pattern_idx == 0; pattern_idx = pattern_len; end
            f0 = this.params.freq_channels(this.hop_pattern(pattern_idx)).f0;
            f1 = this.params.freq_channels(this.hop_pattern(pattern_idx)).f1;

            sq_one = zeros(1, this.params.simstep_sz);
            sq_zero = zeros(1, this.params.simstep_sz);

            stop = length(details.time_slice);
            sq_one((delay_samples + 1):stop) = square(2 * pi * f1 * details.time_slice((delay_samples + 1):stop));
            sq_zero((delay_samples + 1):stop) = square(2 * pi * f0 * details.time_slice((delay_samples + 1):stop));
            
            if (this.curr_step > 1) && (delay_samples > 0)
                prev_pattern_idx = mod(details.prev_steps_symb, pattern_len);
                if prev_pattern_idx == 0; prev_pattern_idx = pattern_len; end
                prev_f0 = this.params.freq_channels(this.hop_pattern(prev_pattern_idx)).f0;
                prev_f1 = this.params.freq_channels(this.hop_pattern(prev_pattern_idx)).f1;
                sq_one(1:delay_samples) = square(2 * pi * prev_f1 * details.time_slice(1:delay_samples));
                sq_zero(1:delay_samples) = square(2 * pi * prev_f0 * details.time_slice(1:delay_samples));
            end
            
            all_ones = sq_one .* details.carrier_slice .* ...
                this.vanatta_gain(this.params.num_elements, ...
                this.params.Fc, this.params.Fc + f1);
            all_zeros = sq_zero .* details.carrier_slice .* ...
                this.vanatta_gain(this.params.num_elements, ...
                this.params.Fc, this.params.Fc + f0);

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
                this bherber.thesis.tags.FreqHopTag
                symb_num double {mustBePositive}
                signal (1, :)
            end

            time_slice = (0:(1 / this.params.Fs):(double(this.params.symb_sz) * (1 / this.params.Fs) - (1 / this.params.Fs))) ...
                + ((symb_num - 1) * double(this.params.symb_sz) * (1 / this.params.Fs));

            carrier_slice = this.params.amplitude * cos(2 * pi * this.params.Fc * time_slice);

            subband_idx = mod(symb_num, length(this.hop_pattern));
            if subband_idx == 0; subband_idx = length(this.hop_pattern); end
            subband = this.freq_subbands(this.hop_pattern(subband_idx));

            sq_one = square(2 * pi * subband.f1 * time_slice);
            sq_zero = square(2 * pi * subband.f0 * time_slice);
            
            all_ones = sq_one .* carrier_slice;
            all_zeros = sq_zero .* carrier_slice;

            % 1. Freq mix
            mixed_one = signal .* all_ones;
            mixed_zero = signal .* all_zeros;

            % 2. Filter
            filtered_one = lowpass(mixed_one, 2 * subband.f1, this.params.Fs);
            filtered_zero = lowpass(mixed_zero, 2 * subband.f0, this.params.Fs);

            % 3. Mean
            correlated_one = trapz(filtered_one);
            correlated_zero = trapz(filtered_zero);

            % 4. Combine Streams
            combined_streams = abs(correlated_one) - abs(correlated_zero);

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
