classdef RepititionChannelCoding
    %REPITION_CHANNEL_CODING unit that implements 2k+1 channel coding.
    %   Detailed explanation goes here

    properties
        repition_rate double {mustBePositive}
    end

    methods
        function obj = RepititionChannelCoding(k)
            %REPITION_CHANNEL_CODING constructor
            %   Implements 2k+1 repitition channel coding scheme.
            obj.repition_rate = k;
        end

        function res = encode(this, bitstring)
            %ENCODE a given bitstring using 2k+1 channel encoding.

            arguments
                this bherber.thesis.channel_coding.RepititionChannelCoding
                bitstring (1, :) {mustBeInRange(bitstring, 0, 1, "inclusive")}
            end

            res = repelem(bitstring, 2 * this.repition_rate + 1);
        end

        function res = decode(this, bitstring)
            %DECODE a given bitstring using 2k+1 channel encoding.
            %   Using Majority vote, decode a given bitstring.

            arguments
                this bherber.thesis.channel_coding.RepititionChannelCoding
                bitstring (1, :) {mustBeInRange(bitstring, 0, 1, "inclusive")}
            end

            vals = reshape(bitstring, [], length(bitstring) / (2 * this.repition_rate + 1));
            res = (mean(vals, 1) >= 0.5);
        end
    end
end