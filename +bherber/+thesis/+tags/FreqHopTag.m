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
        freq_subbands (:, :) double
    end

    %% Public Instance Methods

    methods (Access = public)

        function this = FreqHopTag(x, y, z, mode, bits, sim_params, frequency_subbands, hop_pattern, options)
            %FREQ_HOP_TAG constructor
            %   Create a tag a specific distance from a basestation in 3D
            %   space.
            
            arguments
                x (1, 1) double {mustBeScalarOrEmpty} % x-position of tag relative to basestation
                y (1, 1) double {mustBeScalarOrEmpty} % y-position of tag relative to basestation
                z (1, 1) double {mustBeScalarOrEmpty} % z-position of tag relative to basestation
                mode (1, 1) bherber.thesis.TagType % Modulation mode
                bits (1, :) {mustBeNonmissing, mustBeInRange(bits, 0, 1)} % bitstream
                sim_params (1, 1) bherber.thesis.SimulationConstants % Simulation parameters
                frequency_subbands (:, :) double
                hop_pattern (1, :)
                options.snr_db (1, 1) double = NaN;
                options.complex_noise (1, 1) logical = false;
            end

            this@bherber.thesis.tags.Tag(x, y, z, mode, bits, sim_params, ...
                snr_db=options.snr_db, complex_noise=options.complex_noise);
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
            subband = reshape(this.params.freq_channels(this.hop_pattern(pattern_idx), :), [], 1);

            pathloss = bherber.thesis.PathLoss.amplitude_factor(this.distance, this.params);
            gain = this.vanatta_gain(this.params.num_elements, this.params.Fc, this.params.Fc);

            square_wvs = square(2 .* pi .* subband .* details.time_slice);
            all_symbol_values = square_wvs .* details.carrier_slice .* pathloss .* gain;

            res = zeros(size(details.data));
            for idx = 1:length(res)
                res(idx) = all_symbol_values(details.data(idx), idx);
            end

            if (this.curr_step > 1) && (delay_samples > 0)
                prev_pattern_idx = mod(details.prev_steps_symb, pattern_len);
                if prev_pattern_idx == 0; prev_pattern_idx = pattern_len; end
                prev_subband = reshape(this.params.freq_channels(this.hop_pattern(prev_pattern_idx), :), [], 1);

                prev_square_wvs = square(2 .* pi .* prev_subband .* details.time_slice);
                prev_all_symbol_values = prev_square_wvs .* details.carrier_slice .* pathloss .* gain;
    
                for idx = 1:delay_samples
                    res(idx) = prev_all_symbol_values(details.data(idx), idx);
                end
            end

            if ~isnan(this.snr_db)
                Tsample = 1 / this.params.Fs;
                Tsymbol = 1 / this.params.symb_freq;
                snr_converted = bherber.thesis.channel_models.AWGNChannelModel.EbNo_to_snr(...
                    this.snr_db, Tsample, Tsymbol, this.params.m_ary_modulation, this.complex_noise);
                linear_snr = 10 ^ (snr_converted / 10);
                expected_amplitude = this.params.amplitude * pathloss * gain;
                linear_carrier_power = (expected_amplitude ^ 2) / 4; % ASSUMPTION: Uniform Random Data
                signal_noise = bherber.thesis.channel_models.AWGNChannelModel.noise(...
                    length(details.time_slice), linear_carrier_power, linear_snr, this.complex_noise) + ...
                    bherber.thesis.channel_models.AWGNChannelModel.noise(...
                        length(details.time_slice), linear_carrier_power, linear_snr, this.complex_noise);
                res = res + signal_noise;
            end

            res = res * pathloss;
        end
    
% ----------------------------------------------------------------------- %

        function points = constellation_point(this, symb_num, signal)
            arguments
                this (1, 1) bherber.thesis.tags.FreqHopTag
                symb_num (1, 1) double {mustBePositive}
                signal (1, :)
            end

            time_slice = (0:(1 / this.params.Fs):(double(this.params.symb_sz) * (1 / this.params.Fs) - (1 / this.params.Fs))) ...
                + ((symb_num - 1) * double(this.params.symb_sz) * (1 / this.params.Fs));

            carrier_slice = this.params.amplitude * cos(2 * pi * this.params.Fc * time_slice);

            subband_idx = mod(symb_num, length(this.hop_pattern));
            if subband_idx == 0; subband_idx = length(this.hop_pattern); end
            subband = reshape(this.freq_subbands(this.hop_pattern(subband_idx), :), [], 1);

            square_wvs = square(2 .* pi .* subband .* time_slice);
            all_symbol_values = square_wvs .* carrier_slice;

            % 1. Freq mix
            mixed = signal .* all_symbol_values;

            % 2. Correlate
            correlated = abs(trapz(mixed, 2));
            points = correlated;

        end

% ----------------------------------------------------------------------- %
    
        function res_bits = demodulate(this, symb_num, signal)
            arguments
                this (1, 1) bherber.thesis.tags.FreqHopTag
                symb_num (1, 1) double {mustBePositive}
                signal (1, :)
            end

            correlated = this.constellation_point(symb_num, signal);

            % 5. Decide
            [~, idx] = max(correlated);
            strbits = dec2bin(this.gray2dec(idx - 1), log2(this.params.m_ary_modulation));
            res_bits = this.str2bin(strbits);
        end
    end
end
