classdef AWGNChannelModel < bherber.thesis.channel_models.ChannelModel
    %ADDITIVE_GAUSSIAN_WHITE_NOISE Channel Model.

    properties
        signal
        is_complex;
        snr_db;
        signal_power_db;
        snr_def;
    end

    methods
        function obj = AWGNChannelModel(signal, options)
            arguments
                signal (1, :)
                options.complex (1, 1) = true;
                options.snr_db (1, 1) double = NaN;
                options.signal_power_db (1, 1) double = NaN;
            end

            if isnan(options.signal_power_db)
                assert(~isnan(options.snr_db))
            else
                assert(~isnan(options.signal_power_db))
            end

            obj.signal = signal;
            obj.is_complex = options.complex;
            obj.snr_db = options.snr_db;
            obj.signal_power_db = options.signal_power_db;
        end

        function res = affect_signal(this)
        %AFFECT_SIGNAL with additive Gaussian white noise.

            if isnan(this.signal_power_db)
                assert(~isnan(this.snr_db))
                total_power = sum(this.signal .^2) / numel(this.signal);
            else
                assert(~isnan(this.signal_power_db))
                total_power = 10^(this.signal_power_db / 10);
            end

            linear_snr = 10^(this.snr_db / 10);

            res = this.signal + bherber.thesis.channel_models.AWGNChannelModel.noise(...
                length(this.signal), total_power, linear_snr, this.is_complex);
        end
    end

    methods (Static)
        function vals = noise(n, linear_signal_power, linear_snr, is_complex)
        %AWGN_CHANNEL adding complex gaussian white noise to a signal.
        % Given an snr in db, apply complex white noise to a given signal.
        % Adapted from the built-in toolbox comm.awgn() function.

        arguments
            n (1, 1) double {mustBePositive}
            linear_signal_power (1, 1) double {mustBePositive}
            linear_snr (1, 1) double {mustBePositive}
            is_complex (1, 1) logical
        end
        
            noise = linear_signal_power / linear_snr;
            if is_complex
                vals = (sqrt(noise / 2) * (randn(1, n) + 1j * randn(1, n)));
            else
                vals = (sqrt(noise) * randn(1, n));
            end
        end

        function snr_db = EbNo_to_snr(EbNo, Tsample, Tsymbol, modulation_order, is_complex)
            %EBNO_TO_SNR conversion.
            %   Given a Eb/No value in dB, convert it to an equivalent SNR
            %   value (dB) given the sample and symbol periods, modulation
            %   order, as well as if the noise to be generated should be
            %   complex.

            arguments
                EbNo (1, 1) double
                Tsample (1, 1) double {mustBePositive}
                Tsymbol (1, 1) double {mustBePositive}
                modulation_order (1, 1) double {mustBePositive}
                is_complex (1, 1) logical
            end

            bits_per_symb = log2(modulation_order);
            linear_EbNo = 10 ^ (EbNo / 10);
            if is_complex
                snr_linear = bits_per_symb * linear_EbNo * (Tsample / Tsymbol);
            else
                snr_linear = 0.5 * bits_per_symb * linear_EbNo * (Tsample / Tsymbol);
            end

            snr_db = 10 * log10(snr_linear);

        end
    end
end