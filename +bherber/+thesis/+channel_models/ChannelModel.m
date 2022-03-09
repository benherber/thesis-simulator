classdef ChannelModel
    %CHANNEL_MODEL abstract class.

    methods (Abstract)
        affect_signal(this)
    end

    methods (Abstract, Static)
        noise(n, linear_signal_power, linear_snr, is_complex)
    end
end