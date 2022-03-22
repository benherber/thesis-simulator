% ----------------------------- %
% On-Off Keying Backscatter Tag %
% Name   : OOKTag.m             %
% Author : Benjamin Herber      %
% Date   : Spring 2021          %
% ----------------------------- %

classdef OOKTag < bherber.thesis.tags.Tag
    %OOK_TAG model
    %   Backscatter tag model of a passive Van Atta configuration using
    %   On-Off Keying to communicate

    %% Public Instance Methods

    methods (Access = public)

        function this = OOKTag(x, y, z, mode, bits, sim_params, options)
            %OOK_TAG constructor
            %   Create a tag a specific distance from a basestation in 3D
            %   space.
            
            arguments
                x double {mustBeScalarOrEmpty} % x-position of tag relative to basestation
                y double {mustBeScalarOrEmpty} % y-position of tag relative to basestation
                z double {mustBeScalarOrEmpty} % z-position of tag relative to basestation
                mode bherber.thesis.TagType % Modulation mode
                bits (1, :) {mustBeNonmissing, mustBeInRange(bits, 0, 1)} % bitstream
                sim_params bherber.thesis.SimulationConstants % Simulation parameters
                options.snr_db = NaN;
                options.complex_noise = false;
            end

            this@bherber.thesis.tags.Tag(x, y, z, mode, bits, sim_params, ...
                snr_db=options.snr_db, complex_noise=options.complex_noise);

        end

% ----------------------------------------------------------------------- %
    
        function res = step(this)
            %STEP one simulation step using OOK modulation.
    
            if ~this.is_active
                    res = zeros(1, this.params.simstep_sz);
                    return
            end
    
            details = this.step_data();
            pathloss = bherber.thesis.PathLoss.amplitude_factor(this.distance, this.params);
            gain = this.vanatta_gain(this.params.num_elements, this.params.Fc, this.params.Fc);
            res = details.carrier_slice .* pathloss .* ...
                (this.params.m_ary_amplitudes(details.data) / this.params.amplitude) .* gain;

            if ~isnan(this.snr_db)
                linear_snr = 10 ^ (this.snr_db / 10);
                expected_amplitude = this.params.amplitude * pathloss * gain;
                linear_carrier_power = (expected_amplitude ^ 2) / 4; % ASSUMPTION: Uniform Random Data
                signal_noise = bherber.thesis.channel_models.AWGNChannelModel.noise(...
                    length(details.time_slice), linear_carrier_power, linear_snr, this.complex_noise);
                res = res + signal_noise;
            end

            res = res * pathloss;
        end

% ----------------------------------------------------------------------- %
    
        function point = constellation_point(this, symb_num, signal)
            arguments
                this bherber.thesis.tags.OOKTag
                symb_num double {mustBePositive}
                signal (1, :)
            end

            time_slice = (0:(1 / this.params.Fs):(double(this.params.symb_sz) * (1 / this.params.Fs) - (1 / this.params.Fs))) ...
                + ((symb_num - 1) * double(this.params.symb_sz) * (1 / this.params.Fs));

%             global demod_time
%             demod_time = [demod_time, time_slice];

            carrier_slice = this.params.amplitude * cos(2 * pi * this.params.Fc * time_slice);

            % 1. Freq mix
            mixed_ook = (signal .* carrier_slice);

            % 2. Filter
            filtered = lowpass(mixed_ook, this.params.symb_freq * 2, this.params.Fs);

            % 3. Mean Value
%             correlated_ook = mean(filtered);
%             correlated_ook = rms(filtered) .^ 2;
%             correlated_ook = sum(abs(filtered) .^ 2) / numel(filtered);
            correlated_ook = trapz(filtered);
            point = correlated_ook;

%             gca;
%             hold on
%             ramp = [];
%             for idx = 1:100:(length(time_slice) - 1)
%                 ramp = [ramp, (1 / this.params.Fs) * trapz(filtered(1:idx))];
%             end
%             if isnan(this.snr_db)
%                 plot(gca, ramp, "-b");
%             else
%                 plot(gca, ramp, "-r");
%             end
        end
    
% ----------------------------------------------------------------------- %
    
        function res_bits = demodulate(this, symb_num, signal)
            arguments
                this bherber.thesis.tags.OOKTag
                symb_num double {mustBePositive}
                signal (1, :)
            end

            correlated_ook = this.constellation_point(symb_num, signal);
            
            % 4. Decide
            expected_amplitudes = this.params.m_ary_amplitudes .* ...
                ((bherber.thesis.PathLoss.amplitude_factor(this.distance, this.params) .^ 2)) .* ...
                 this.vanatta_gain(this.params.num_elements, this.params.Fc, this.params.Fc);
            filtered_amp = (expected_amplitudes .* this.params.amplitude) / 2;
            lambda = double(this.params.symb_sz - 1) .* filtered_amp;
            [~, idx] = min(abs(correlated_ook - lambda));
            strbits = dec2bin(this.gray2dec(idx - 1), log2(this.params.m_ary_modulation));
            res_bits = this.str2bin(strbits);

%             global ook_mixed; global ook_filtered; global ook_correlated;
%             ook_mixed = [ook_mixed, mixed_ook];
%             ook_filtered = [ook_filtered, filtered];
%             ook_correlated = [ook_correlated, repelem(correlated_ook, this.params.symb_sz)];
        end
    end
end
