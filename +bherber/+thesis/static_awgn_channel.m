function res = static_awgn_channel(signal, params, options)
%STATIC_AWGN_CHANNEL adding complex gaussian white noise to a signal.
%  Given an snr in db, apply complex white noise to a given signal which
%  is statically calculated based off the amplitude of the signal and
%  modulation type.

    arguments
        signal (1, :)
        params bherber.thesis.SimulationConstants
        options.mode bherber.thesis.TagType = bherber.thesis.TagType.NULL;
        options.distance double = NaN;
        options.snr_db double = NaN;
        options.noise_power_db double = NaN;
    end
    
    if options.mode ~= bherber.thesis.TagType.NULL
        assert(~isnan(options.snr_db));
        assert(~isnan(options.distance));
        switch options.mode
            case bherber.thesis.TagType.OOK
                total_power = (params.amplitude .^2) / 2;
            case {  bherber.thesis.TagType.FSK_LO, ...
                    bherber.thesis.TagType.FSK_HI, ...
                    bherber.thesis.TagType.FREQ_HOP  }
                total_power = params.amplitude .^2;
            otherwise
                error("Error -- Unsupported Tag Type: %s", string(options.mode));
        end
        if options.distance > 0
            total_power = total_power * ((params.wavelen / (4 * pi * options.distance)) .^ 2);
        end
        linear_snr = 10 ^ (options.snr_db / 10);
        variance = sqrt(total_power / linear_snr);
        res = signal + sqrt(variance / 2) * (randn(size(signal)) + 1j * randn(size(signal)));
    else
        assert(~isnan(options.noise_power_db))
        noise = 10 ^ (options.noise_power_db / 10);
        res = signal + (sqrt(noise / 2) * (randn(size(signal)) + 1j * randn(size(signal))));
    end
    
end