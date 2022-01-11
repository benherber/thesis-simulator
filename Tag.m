% ------------------------ %
% Back Scatter Tag         %
% Name   : Tag.m           %
% Author : Benjamin Herber %
% Date   : Fall 2021       %
% ------------------------ %

classdef Tag < handle
    %TAG model
    %   Backscatter tag model of a passive Van Atta configuration

    %%  Private Properties

    properties (GetAccess = public, SetAccess = private)
        x double {mustBeScalarOrEmpty} % x-position of tag relative to basestation
        y double {mustBeScalarOrEmpty} % y-position of tag relative to basestation
        z double {mustBeScalarOrEmpty} % z-position of tag relative to basestation
        mode TagType % Modulation mode
        curr_step int64 {mustBeScalarOrEmpty} % Current time-slice in modulation
        carrier (:, 1) {mustBeFinite, mustBeNonmissing} % Delayed carrier signal
        data (:, 1) {mustBeNonmissing, mustBeInRange(data, 0, 1)} % Delayed data signal
        time (:, 1) {mustBeReal, mustBeFinite, mustBeNonmissing} % Time signal
        channel function_handle % Propagation channel
        params SimulationConstants % Simulation parameters
    end

    %%  Dependent Properties

    properties (Dependent)
        distance % Tag distance from basestation
    end

    %%  Private Static Methods
    methods (Static, Access = private)
        function phase = find_phase(theta, wavelen, spacing, n)
            % FIND_PHASE Find phase shift of an element in a phased array.
            %   Given an impinging angle, wavelength, inter-element
            %   spacing, and element number in a phased array, calculate
            %   the phase shift.

            phase = ((2.0 * pi) / wavelen) * spacing * n * cos(theta);

        end

        % ----------------------------------------------------------------------- %

        function res = vanatta_gain(num_elements, fc, fb)
            %VANATTA_GAIN of a backscatter tag in a Van Atta Configuration.
            %   Calculate gain of a backscatter tag in a Van Atta
            %   Configuration for given input and output frequencies.

            tmp = 0;
            C = physconst("Lightspeed");
            lambda_c = C / fc;
            lambda_b = C / fb;
            spacing = lambda_c / 2;

            for idx = 1:num_elements
                tmp = tmp + exp(1i * ...
                    (Tag.find_phase(0, lambda_b, spacing, (idx - 1)) ...
                    - Tag.find_phase(0, lambda_c, spacing, (idx - 1))));
            end

            res = tmp;
        end

        % ----------------------------------------------------------------------- %

        function res = friis_path_loss(signal, distance, sim_params)
            %FRIIS_PATH_LOSS of signal through free space
            %   Using Friis Formula for path loss through free space, calculate
            %   the remaining wave.

            if distance == 0
                res = signal;
            else
                % Calculate power loss
                power_loss = sim_params.wavelen / ...
                    (16 * pi * pi * distance * distance);

                % Calculate the resultant wave
                res = signal * sqrt(2) * sqrt(power_loss);
            end
        end

    end

    %%  Private Methods
    methods (Access = private)

        function res = delay(this, signal)
            %DELAY Delay a signal
            %   Delay a signal based on the tag's distance from the
            %   basestation.

            time_delay = this.distance / physconst("Lightspeed");
            num_samples = time_delay * this.params.Fs;
            res = [zeros(int32(num_samples), 1).', signal];
        end

        % ----------------------------------------------------------------------- %

        function [f1, f0] = freq(this)
            %FREQ determine current operating frequency

            if this.mode == TagType.FSK_LO
                f1 = this.params.fsk_channel0.f1;
                f0 = this.params.fsk_channel0.f0;

            elseif this.mode == TagType.FSK_HI
                f1 = this.params.fsk_channel1.f1;
                f0 = this.params.fsk_channel1.f0;

            elseif this.mode == TagType.FSK_RANDOM

                if logical(randi([0, 1]))
                    f1 = this.params.fsk_channel0.f1;
                    f0 = this.params.fsk_channel0.f0;
                else
                    f1 = this.params.fsk_channel1.f1;
                    f0 = this.params.fsk_channel1.f0;
                end

            else
                error("UNSUPPORTED MODE");
            end

        end
    end

    %% Public Methods

    methods

        function this = Tag(x, y, z, mode, time, carrier, data, channel, sim_params)
            %TAG constructor
            %   Create a tag a specific distance from a basestation in 3D
            %   space.

            this.params = sim_params;
            this.x = x;
            this.y = y;
            this.z = z;
            this.mode = mode;
            this.curr_step = int64(0);
            this.carrier = this.delay(carrier);
            this.data = this.delay(data);
            this.time = time;
            this.channel = channel;

        end

        % ----------------------------------------------------------------------- %

        function res = get.distance(this)
            %DISTANCE from basestation.
            %   Calculate distance from basestation if it is considered to
            %   be at the origin in 3D space.

            res = sqrt((this.x)^2 + (this.y)^2 + (this.z)^2);
        end

        function res = step(this, n)
            %STEP the object 'n' simulation frame(s) forward

            arguments
                this Tag
                n {mustBeFinite, mustBeScalarOrEmpty, mustBePositive} = 1
            end

            start_pt = uint64(this.curr_step * this.params.simstep_sz + 1);
            stop_pt = uint64(start_pt + (n * this.params.simstep_sz) - 1);

            if stop_pt > length(this.time)
                error("SIMULATION RAN OUT OF BOUNDS");
            end

            cut_time = this.time(start_pt:stop_pt);
            cut_carrier = this.carrier(start_pt:stop_pt);
            cut_data = this.data(start_pt:stop_pt);

            cut_carrier = Tag.friis_path_loss(this.channel(cut_carrier), ...
                this.distance, this.params);
            if this.mode == TagType.OOK
                res = Tag.modulate_by_ook(cut_carrier, cut_data, this.params);
            else
                [f1, f0] = this.freq();
                res = Tag.modulate_by_fsk(cut_time, cut_carrier, ...
                    cut_data, this.params, f1, f0);
            end
            res = Tag.friis_path_loss(this.channel(res), ...
                this.distance, this.params);

            this.curr_step = this.curr_step + n;
        end

    end

    %% Static Public Methods
    methods (Static)

        function ook_modulation = modulate_by_ook(carrier, data, params)
            % MODULATE_BY_OOK Modulate a given signal by on-off keying.
            %   Given carrier and data signals, modulate the carrier wave
            %   according to on-off keying as a backscatter tag would
            %   through a given channel.

            % 2. Modulate
            ook_modulation = carrier .* data .* ...
                Tag.vanatta_gain(params.num_elements, params.Fc, params.Fc);

        end

        % ----------------------------------------------------------------------- %

        function f_modulation = modulate_by_fsk(t, carrier, data, params, f1, f0)
            % MODULATE_BY_FSK a given signal.
            %   Given carrier and data signals, modulate the carrier wave
            %   according to frequency-shift keying as a backscatter tag
            %   would through a given channel.

            % Frequency modulation
            f_modulation = zeros(size(data));
            sq_one = square((f1) * 2 * pi * t);
            sq_zero = square((f0) * 2 * pi * t);
            all_ones = sq_one .* carrier * ...
                Tag.vanatta_gain(params.num_elements, params.Fc, params.Fc + f1);
            all_zeros = sq_zero .* carrier * ...
                Tag.vanatta_gain(params.num_elements, params.Fc, params.Fc + f0);

            for idx = 1:length(f_modulation)

                if data(idx) == 1
                    f_modulation(idx) = all_ones(idx);
                else
                    f_modulation(idx) = all_zeros(idx);
                end

            end

        end
        
        % ----------------------------------------------------------------------- %

        function res_bits = ook_demodulate(signal, carrier, time)
            %OOK_DEMODULATE a signal
            %   Demodulate a signal modulated by on-off keying.

            arguments
                signal (:, 1) {mustBeFinite, mustBeNonmissing}
                carrier (:, 1) {mustBeFinite, mustBeNonmissing}
                time (:, 1) {mustBeFinite, mustBeReal}
            end

            % 1. Freq mix
            mixed_ook = signal .* carrier;

            % 2. Integrate
            correlated_ook = trapz(time, mixed_ook);

            % 3. Decide
            lambda = trapz(time, carrier) / 2;
            if correlated_ook > lambda
                res_bits = 1;
            else
                res_bits = 0;
            end
        end

        function res_bits = fsk_demodulate(signal, carrier, time, f1, f0)
            %FSK_DEMODULATE a signal
            %   Demodulate a signal modulated by frequency-shift keying with
            %   frequencies f1 and f0 corresponding to 'on' & 'off'
            %   respectively.

            arguments
                signal (:, 1) {mustBeFinite, mustBeNonmissing}
                carrier (:, 1) {mustBeFinite, mustBeNonmissing}
                time (:, 1) {mustBeFinite, mustBeReal}
                f1 {mustBeNumeric, mustBeReal}
                f0 {mustBeNumeric, mustBeReal}
            end

            % 1. Freq mix
            all_ones = square((f1) * 2 * pi * time) .* carrier;
            all_zeros = square((f0) * 2 * pi * time) .* carrier;
            mixed_one = signal .* all_ones;
            mixed_zero = signal .* all_zeros;

            % 2. Integrate and sharpen
            correlated_one = trapz(time, mixed_one);
            correlated_zero = trapz(time, mixed_zero);

            % 3. Combine Streams
            combined_streams = abs(correlated_one) - abs(correlated_zero);

            % 4. Decide
            lambda = 0;
            if combined_streams > lambda
                res_bits = 1;
            else
                res_bits = 0;
            end

        end

    end

end
