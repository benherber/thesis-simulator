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
        mode bherber.thesis.TagType % Modulation mode
        curr_step int64 {mustBeScalarOrEmpty} % Current time-slice in modulation
        carrier (:, 1) {mustBeFinite, mustBeNonmissing} % Delayed carrier signal
        data (:, 1) {mustBeNonmissing, mustBeInRange(data, 0, 1)} % Delayed data signal
        time (:, 1) {mustBeReal, mustBeFinite, mustBeNonmissing} % Time signal
        channel function_handle % Propagation channel
        fh_pattern (:, 1) = [];
        is_active logical = true;
        square_one
        square_zero
        params bherber.thesis.SimulationConstants % Simulation parameters
    end

    %%  Dependent Properties

    properties (Dependent)
        distance % Tag distance from basestation
    end

    %%  Private Static Methods
    methods (Static)
        function phase = find_phase(theta, wavelen, spacing, n)
            % FIND_PHASE Find phase shift of an element in a phased array.
            %   Given an impinging angle, wavelength, inter-element
            %   spacing, and element number in a phased array, calculate
            %   the phase shift.

            phase = ((2.0 * pi) / wavelen) * spacing * n * cos(theta);  % FIXME: look at again

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
                    (bherber.thesis.Tag.find_phase(0, lambda_b, spacing, (idx - 1)) ...
                    - bherber.thesis.Tag.find_phase(0, lambda_c, spacing, (idx - 1))));
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
                power_loss = sim_params.wavelen / (4 * pi * distance);

                % Calculate the resultant wave
                res = signal * power_loss;
            end
        end

    end

    %%  Private Methods
    methods (Access = private)

        function res = delay(this, signal)
            %DELAY Delay a signal
            %   Delay a signal based on the tag's distance from the
            %   basestation.

            arguments
                this bherber.thesis.Tag
                signal (:, 1)
            end

            time_delay = this.distance / physconst("Lightspeed");
            num_samples = time_delay * this.params.Fs;
            res = [zeros(int32(num_samples), 1); signal];
        end
    end

    %% Public Methods

    methods

        function this = Tag(x, y, z, mode, time, carrier, data, channel, sim_params)
            %TAG constructor
            %   Create a tag a specific distance from a basestation in 3D
            %   space.
            
            arguments
                x double {mustBeScalarOrEmpty} % x-position of tag relative to basestation
                y double {mustBeScalarOrEmpty} % y-position of tag relative to basestation
                z double {mustBeScalarOrEmpty} % z-position of tag relative to basestation
                mode bherber.thesis.TagType % Modulation mode
                time (:, 1) {mustBeReal, mustBeFinite, mustBeNonmissing} % Time signal
                carrier (:, 1) {mustBeFinite, mustBeNonmissing} % Delayed carrier signal
                data (:, 1) {mustBeNonmissing, mustBeInRange(data, 0, 1)} % Delayed data signal
                channel function_handle % Propagation channel
                sim_params bherber.thesis.SimulationConstants % Simulation parameters
            end

            this.params = sim_params;
            this.x = x;
            this.y = y;
            this.z = z;
            this.mode = mode;
            this.curr_step = int64(0);
            this.carrier = carrier; % NO DELAY
            this.data = data; % NO DELAY
            this.time = time;
            this.channel = channel;

        end

% ----------------------------------------------------------------------- %
        
        function res = generate_hoppattern(this)
            if this.mode == bherber.thesis.TagType.FREQ_HOP
                res = randi([1, this.params.num_channels], 1, ...
                    this.params.pattern_len);
                this.fh_pattern = res;
            else
                error("Tag is not in FREQ_HOP mode.");
            end
        end

% ----------------------------------------------------------------------- %

        function activate(this)
            this.is_active = true;
        end

        function deactivate(this)
            this.is_active = false;
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
                this bherber.thesis.Tag
                n {mustBeFinite, mustBeScalarOrEmpty, mustBePositive} = 1
            end

            if ~this.is_active
                res = zeros(1, this.params.simstep_sz);
                return
            end

            if this.curr_step == 0
                if this.mode ~= bherber.thesis.TagType.OOK
                    if this.mode == bherber.thesis.TagType.FREQ_HOP
                        shaped_time = reshape(this.time, [round(length(this.time) / (this.params.num_symbs + 1)), ...
                            this.params.num_symbs + 1]);
                        sq_one = zeros(round(length(this.time) / (this.params.num_symbs + 1)), ...
                            this.params.num_symbs + 1);
                        sq_zero = zeros(round(length(this.time) / (this.params.num_symbs + 1)), ...
                            this.params.num_symbs + 1);
                        for idx = 1:this.params.num_symbs
                            channel_idx = mod(idx, this.params.pattern_len);
                            if channel_idx == 0; channel_idx = this.params.pattern_len; end
                            f0 = this.params.freq_channels(this.fh_pattern(channel_idx)).f0;
                            f1 = this.params.freq_channels(this.fh_pattern(channel_idx)).f1;
                            sq_one(:, idx) = square(f1 * 2 * pi * shaped_time(:, idx));
                            sq_zero(:, idx) = square(f0 * 2 * pi * shaped_time(:, idx));
                        end
                        this.square_one = sq_one(:); % this.delay(sq_one(:));
                        this.square_zero = sq_zero(:); % this.delay(sq_zero(:));
                    else
                        if this.mode == bherber.thesis.TagType.FSK_LO
                            f0 = this.params.freq_channels(1).f0;
                            f1 = this.params.freq_channels(1).f1;
                        elseif this.mode == bherber.thesis.TagType.FSK_HI
                            f0 = this.params.freq_channels(2).f0;
                            f1 = this.params.freq_channels(2).f1;
                        end
                        this.square_one = square(f1 * 2 * pi * this.time); % this.delay(square(f1 * 2 * pi * this.time));
                        this.square_zero = square(f0 * 2 * pi * this.time); % this.delay(square(f0 * 2 * pi * this.time));
                    end
                end
            end

            start_pt = uint64(this.curr_step * this.params.simstep_sz + 1);
            stop_pt = uint64(start_pt + (n * this.params.simstep_sz) - 1);

            if stop_pt > length(this.time)
                error("SIMULATION RAN OUT OF BOUNDS");
            end

            cut_carrier = this.carrier(start_pt:stop_pt);
            cut_data = this.data(start_pt:stop_pt);

%             noisy_carr = this.channel(cut_carrier);
            cut_carrier = bherber.thesis.Tag.friis_path_loss(cut_carrier, ...
            this.distance, this.params);
            if this.mode == bherber.thesis.TagType.OOK
                res = bherber.thesis.Tag.modulate_by_ook(cut_carrier, cut_data, this.params);
            else
                cut_sq_one = this.square_one(start_pt:stop_pt);
                cut_sq_zero = this.square_zero(start_pt:stop_pt);
                if this.mode == bherber.thesis.TagType.FREQ_HOP
                    symb_number = ceil(this.curr_step / this.params.sim_sym_ratio);
                    pattern_idx = mod(symb_number, this.params.pattern_len);
                    if pattern_idx == 0; pattern_idx = this.params.pattern_len; end
                    f0 = this.params.freq_channels(this.fh_pattern(pattern_idx)).f0;
                    f1 = this.params.freq_channels(this.fh_pattern(pattern_idx)).f1;
                else
                    if this.mode == bherber.thesis.TagType.FSK_LO
                        channel_idx = 1;
                    else
                        channel_idx = 2;
                    end
                    f0 = this.params.freq_channels(channel_idx).f0;
                    f1 = this.params.freq_channels(channel_idx).f1;
                end
                res = bherber.thesis.Tag.modulate_by_fsk(cut_carrier, cut_data, ...
                    this.params, cut_sq_one, cut_sq_zero, f1, f0);
            end
%             noisy_res = this.channel(res);
            res = bherber.thesis.Tag.friis_path_loss(res, ...
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

            arguments
                carrier (1, :)
                data (1, :)
                params bherber.thesis.SimulationConstants
            end

            ook_modulation = carrier .* data .* ...
                bherber.thesis.Tag.vanatta_gain(params.num_elements, params.Fc, params.Fc);

        end

        % ----------------------------------------------------------------------- %

        function f_modulation = modulate_by_fsk(carrier, data, params, sq_one, sq_zero, f1, f0)
            % MODULATE_BY_FSK a given signal.
            %   Given carrier and data signals, modulate the carrier wave
            %   according to frequency-shift keying as a backscatter tag
            %   would through a given channel.

            arguments
                carrier (1, :)
                data (1, :)
                params bherber.thesis.SimulationConstants
                sq_one (1, :)
                sq_zero (1, :)
                f1 {mustBePositive}
                f0 {mustBePositive}
            end

            % Frequency modulation
            f_modulation = zeros(size(data));
            all_ones = sq_one .* carrier * ...
                bherber.thesis.Tag.vanatta_gain(params.num_elements, params.Fc, params.Fc + f1);
            all_zeros = sq_zero .* carrier * ...
                bherber.thesis.Tag.vanatta_gain(params.num_elements, params.Fc, params.Fc + f0);

            for idx = 1:length(f_modulation)
                if data(idx) == 1
                    f_modulation(idx) = all_ones(idx);
                else
                    f_modulation(idx) = all_zeros(idx);
                end

            end

        end
        
        % ----------------------------------------------------------------------- %

        function res_bits = ook_demodulate(signal, carrier, time, tag)
            %OOK_DEMODULATE a signal
            %   Demodulate a signal modulated by on-off keying.

            arguments
                signal (1, :)
                carrier (1, :) {mustBeFinite, mustBeNonmissing}
                time (1, :) {mustBeFinite, mustBeReal}
                tag bherber.thesis.Tag
            end

            % 1. Freq mix
            global mixed;
            global sig;
            global filter;
            sig = [sig, signal];
            mixed_ook = (signal .* carrier);
            mixed = [mixed, mixed_ook];

            % 2. Filter
            filtered = lowpass(mixed_ook, tag.params.symb_freq * 2, tag.params.Fs);
            filter = [filter, filtered];

            % 3. Integrate
            correlated_ook = mean(filtered);
%             fprintf("%0.10f\n", abs(correlated_ook));
            global points;
            points = [points; [real(correlated_ook), imag(correlated_ook)]];
            
            % 4. Decide
            expected_amplitude = tag.params.amplitude * ...
                ((tag.params.wavelen / (4 * pi * tag.distance)) ^ 2) * ...
                 tag.vanatta_gain(tag.params.num_elements, tag.params.Fc, tag.params.Fc);
            lambda = (expected_amplitude * tag.params.amplitude) / 4;
            if real(correlated_ook) > lambda
                res_bits = 1;
            else
                res_bits = 0;
            end
        end

        function res_bits = fsk_demodulate(signal, carrier, time, f1, f0, params)
            %FSK_DEMODULATE a signal
            %   Demodulate a signal modulated by frequency-shift keying with
            %   frequencies f1 and f0 corresponding to 'on' & 'off'
            %   respectively.

            arguments
                signal (1, :)
                carrier (1, :) {mustBeFinite, mustBeNonmissing}
                time (1, :) {mustBeFinite, mustBeReal}
                f1 {mustBeNumeric, mustBeReal}
                f0 {mustBeNumeric, mustBeReal}
                params bherber.thesis.SimulationConstants
            end

            % 1. Freq mix
            all_ones = square(f1 * 2 * pi * time) .* carrier;
            all_zeros = square(f0 * 2 * pi * time) .* carrier;
            mixed_one = signal .* all_ones;
            mixed_zero = signal .* all_zeros;

            % 2. Filter
            filtered_one = lowpass(mixed_one, params.Fc / 4, params.Fs);
            filtered_zero = lowpass(mixed_zero, params.Fc / 4, params.Fs);

            % 3. Mean
            correlated_one = mean(filtered_one);
            correlated_zero = mean(filtered_zero);

            % 4. Combine Streams
            combined_streams = correlated_one - correlated_zero;
%             fprintf("%0.2f %0.2fj\n", real(combined_streams), imag(combined_streams));
%             hold on
%             scatter(gca, real(combined_streams), imag(combined_streams), "filled"); 

            % 4. Decide
            lambda = 0;
            if real(combined_streams) > lambda
                res_bits = 1;
            else
                res_bits = 0;
            end

        end

        % --------------------------------------------------------------- %
        function res_bits = fh_demodulate(signal, carrier, time, params)
            %FSK_DEMODULATE a signal
            %   Demodulate a signal modulated by frequency-shift keying with
            %   frequencies f1 and f0 corresponding to 'on' & 'off'
            %   respectively.

            arguments
                signal (1, :)
                carrier (1, :) {mustBeFinite, mustBeNonmissing}
                time (1, :) {mustBeFinite, mustBeReal}
                params bherber.thesis.SimulationConstants
            end

            % 1. DEHOP
            % 2. FSK DEMODULATOR

            % 1. Freq Presence Test
            filter_channels = zeros(params.num_channels, length(signal));
            for idx = 1:params.num_channels
                lo = params.freq_channels(idx).f0 - (params.interchannel_spacing / 2);
                hi = params.freq_channels(idx).f1 + (params.interchannel_spacing / 2);
                filter_channels(idx, :) = bandpass(signal, [lo, hi], params.Fs);
            end
            
            res = NaN(1, params.num_channels);
            all_ones = square((f1) * 2 * pi * time) .* carrier;
            all_zeros = square((f0) * 2 * pi * time) .* carrier;
            mixed_one = signal .* all_ones;
            mixed_zero = signal .* all_zeros;

            % 2. Integrate and sharpen
            correlated_one = trapz(time, mixed_one);
            correlated_zero = trapz(time, mixed_zero);

            % 3. Combine Streams
            combined_streams = correlated_one - correlated_zero;
%             fprintf("%0.2f %0.2fj\n", real(combined_streams), imag(combined_streams));

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
