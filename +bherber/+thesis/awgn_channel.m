function res = awgn_channel(signal, options)
%AWGN_CHANNEL adding complex gaussian white noise to a signal.
% Given an snr in db, apply complex white noise to a given signal.
% Adapted from the built-in toolbox comm.awgn() function.

    arguments
        signal (1, :)
        options.complex = true;
        options.snr_db double = NaN;
        options.noise_power_db double = NaN;
    end

    if isnan(options.noise_power_db)
        assert(~isnan(options.snr_db))
        total_power = sum(signal.^2) / numel(signal);
        linear_snr = 10^(options.snr_db / 10);
        noise = total_power / linear_snr; % \sigma^2
    else
        assert(~isnan(options.noise_power_db))
        noise = 10^(options.noise_power_db / 10);
    end

    if options.complex
        noise = (sqrt(noise / 2) * (randn(size(signal)) + 1j * randn(size(signal))));
    else
        noise = signal + (sqrt(noise) * randn(size(signal)));
    end
    res = signal + noise;
end