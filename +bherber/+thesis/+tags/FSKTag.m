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
        freqs (:, 1) {mustBePositive}
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
                frequency_subband (1, :)
                options.snr_db = NaN;
                options.complex_noise = false;
            end

            this@bherber.thesis.tags.Tag(x, y, z, mode, bits, sim_params, ...
                snr_db=options.snr_db, complex_noise=options.complex_noise);
            this.freqs = frequency_subband;

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
            square_wvs = square(2 .* pi .* this.freqs .* details.time_slice);
            
            all_symbol_values = square_wvs .* details.carrier_slice .* gain;

            res = zeros(size(details.data));
            for idx = 1:length(res)
                res(idx) = all_symbol_values(details.data(idx), idx);
            end

            if ~isnan(this.snr_db)
                Tsample = 1 / this.params.Fs;
                Tsymbol = 1/ this.params.symb_freq;
                snr_converted = bherber.thesis.channel_models.AWGNChannelModel.EbNo_to_snr(...
                    this.snr_db, Tsample, Tsymbol, this.params.m_ary_modulation, this.complex_noise);
                linear_snr = 10 ^ (snr_converted / 10);
                expected_amplitude = this.params.amplitude * pathloss * gain;
                linear_carrier_power = (expected_amplitude ^ 2) / 2; % ASSUMPTION: Uniform Random Data
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
                this bherber.thesis.tags.FSKTag
                symb_num double {mustBePositive}
                signal (1, :)
            end

            time_slice = (0:(1 / this.params.Fs):(double(this.params.symb_sz) * (1 / this.params.Fs) - (1 / this.params.Fs))) ...
                + ((symb_num - 1) * double(this.params.symb_sz) * (1 / this.params.Fs));

            carrier_slice = this.params.amplitude * cos(2 * pi * this.params.Fc * time_slice);
            
            square_wvs = square(2 .* pi .* this.freqs .* time_slice);
            all_symbol_values = square_wvs .* carrier_slice;

            % 1. Freq mix
            mixed = signal .* all_symbol_values;

            % 2. Correlate
            correlated = trapz(mixed, 2);
            points = correlated;

        end
    
% ----------------------------------------------------------------------- %
    
    function res_bits = demodulate(this, symb_num, signal)
            arguments
                this bherber.thesis.tags.FSKTag
                symb_num double {mustBePositive}
                signal (1, :)
            end

            correlated = this.constellation_point(symb_num, signal);

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
            [~, idx] = max(correlated);
            strbits = dec2bin(this.gray2dec(idx - 1), log2(this.params.m_ary_modulation));
            res_bits = this.str2bin(strbits);
        end
    end
end
