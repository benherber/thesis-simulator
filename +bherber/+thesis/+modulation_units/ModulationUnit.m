classdef ModulationUnit
    %MODULATION_UNIT Abstract class for any recognized modulation unit.

    properties (Access = protected)
        signal (1, :)
        data (1, :)
        carrier (1, :) {mustBeFinite, mustBeNonmissing}
        params bherber.thesis.SimulationConstants
        distance double
        square_one
        square_zero
        freq_channel struct
    end

    methods
        function obj = ModulationUnit(signal, data, carrier, params, distance, ...
                square_one, square_zero, freq_channel)

            obj.signal = signal;
            obj.data = data;
            obj.carrier = carrier;
            obj.params = params;
            obj.distance = distance;
            obj.square_one = square_one;
            obj.square_zero = square_zero;
            obj.freq_channel = freq_channel;
        end
    end

    methods (Static)
        function res = factory(type, signal, data, carrier, params, distance, square_one, ...
                        square_zero, freq_channels)
            switch type
                case bherber.thesis.TagType.OOK
                    res = bherber.thesis.modulation_units.OOKModulationUnit(...
                        signal, data, carrier, params, distance, [], [], []);
                case bherber.thesis.TagType.FSK_LO
                    res = bherber.thesis.modulation_units.FSKModulationUnit(...
                        signal, data, carrier, params, distance, square_one, ...
                        square_zero, freq_channels(1);
                case bherber.thesis.TagType.FSK_HI
                    res = bherber.thesis.modulation_units.FSKModulationUnit(...
                        signal, data, carrier, params, distance, square_one, ...
                        square_zero, freq_channels(2);
            end
        end
    end

    methods (Abstract)
        modulate(this)
        demodulate(this)
    end
end