function res = awgn_channel(signal, options)
%AWGN_CHANNEL adding complex gaussian white noise to a signal.
% Given an snr in db, apply complex white noise to a given signal.
% Adapted from the built-in toolbox comm.awgn() function.

    arguments
        signal (1, :)
        options.snr_db double = NaN;
        options.noise_power_db double = NaN;
    end

    if isnan(options.noise_power_db)
        assert(~isnan(options.snr_db))
        total_power = sum(abs(signal(:)).^2) / numel(signal(:));
        linear_snr = 10^(options.snr_db / 10);
        noise = total_power / linear_snr;
    else
        assert(~isnan(options.noise_power_db))
        noise = options.noise_power_db;
    end

    res = signal + (sqrt(noise / 2) * (randn(size(signal)) + 1j * randn(size(signal))));
end