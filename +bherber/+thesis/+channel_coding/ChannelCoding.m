classdef ChannelCoding
    %CHANNEL_CODING abstract class.

    methods (Abstract)
        encode(this, bitstring)
        decode(this, bitstring)
    end
end