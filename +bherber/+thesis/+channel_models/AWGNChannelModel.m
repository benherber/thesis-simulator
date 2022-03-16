classdef AWGNChannelModel < bherber.thesis.channel_models.ChannelModel
    %ADDITIVE_GAUSSIAN_WHITE_NOISE Channel Model.

    properties
        signal
        is_complex;
        snr_db;
        signal_power_db;
    end

    methods
        function obj = AWGNChannelModel(signal, options)
            arguments
                signal (1, :)
                options.complex = true;
                options.snr_db double = NaN;
                options.signal_power_db double = NaN;
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
        
            noise = linear_signal_power / linear_snr;
            if is_complex
                vals = (sqrt(noise / 2) * (randn(1, n) + 1j * randn(1, n)));
            else
                vals = (sqrt(noise) * randn(1, n));
            end
        end
    end
end