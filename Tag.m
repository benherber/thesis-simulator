% ------------------------ %
% Back Scatter Tag         %
% Name   : Tag.m           %
% Author : Benjamin Herber %
% Date   : Fall 2021       %
% ------------------------ %

classdef Tag < handle
    %TAG model
    %   Backscatter tag model of a passive Van Atta configuration

    properties
        x % x-position of tag relative to basestation
        y % y-position of tag relative to basestation
        z % z-position of tag relative to basestation
        mode % Frequency mode ('lo', 'hi', or 'random' bands)
    end

    properties (Access = private)
        curr_step % Current time-slice in modulation
        carrier % Delayed carrier signal
        data % Delayed data signal
        time % Time signal
        channel; % Propagation channel
    end

    properties (Dependent)
        distance % Tag distance from basestation
    end

    methods (Access = private)

        function res = delay(this, signal)
            %DELAY Delay a signal
            %   Delay a signal based on the tag's distance from the
            %   basestation.

            time_delay = this.distance / physconst("Lightspeed");
            num_samples = int32(time_delay * Fs);
            padded = [signal; zeros(num_elements, 1)];
            res = delayseq(padded, num_samples);
        end

        function [f1, f0] = freq(this)
            %FREQ determine current operating frequency
            switch this.mode
                case "lo"
                    f1 = FB1;
                    f0 = FB0;

                case "hi"
                    f1 = FB11;
                    f0 = FB01;

                case "random"

                    if logical(randi([0, 1]))
                        f1 = FB1;
                        f0 = FB0;
                    else
                        f1 = FB11;
                        f0 = FB01;
                    end

                otherwise
                    error("UNSUPPORTED MODE");
            end

        end

        function modulated = modulate_by_fsk(this, t, carrier, data, channel)
            % MODULATE_BY_FSK a given signal.
            %   Given carrier and data signals, modulate the carrier wave according to
            %   frequency-shift keying as a backscatter tag would through a given channel.

            % Load Contants

            simulation_constants

            % Frequency modulation

            % 1. Add channel noise headed to tag
            noisy_carrier = friis_path_loss(channel(carrier), this.distance());

            % 2. Modulate
            f_modulation = zeros(size(data));
            sq_one = square((FB1) * 2 * pi * t);
            sq_zero = square((FB0) * 2 * pi * t);
            [f1, f0] = this.freq();
            all_ones = sq_one .* vanatta_gain(noisy_carrier, NUM_ELEMENTS, FC, f1);
            all_zeros = sq_zero .* vanatta_gain(noisy_carrier, NUM_ELEMENTS, FC, f0);

            for idx = 1:length(f_modulation)

                if data(idx) == 1
                    f_modulation(idx) = all_ones(idx);
                else
                    f_modulation(idx) = all_zeros(idx);
                end

            end

            % 3. Add channel noise headed back to basestation
            modulated = friis_path_loss(channel(f_modulation), this.distance());
        end

    end

    methods

        function this = Tag(x, y, z, mode, time, carrier, data, channel)
            %TAG constructor
            %   Create a tag a specific distance from a basestation in 3D
            %   space.

            this.x = x;
            this.y = y;
            this.z = z;
            this.mode = mode;
            this.curr_step = 0;
            this.carrier = this.delay(carrier);
            this.data = this.delay(data);
            this.time = time;
            this.channel = channel;
        end

        function res = get.distance(this)
            %DISTANCE from basestation.
            %   Calculate distance from basestation if it is considered to
            %   be at the origin in 3D space.

            res = sqrt((this.x).^2 + (this.y).^2 + (this.z).^2);
        end

        function [this, res] = step(this, n)
            %STEP the object 'n' simulation frame(s) forward

            if nargin < 2
                n = 1;
            end

            start_pt = this.curr_step * SIMSTEP_SIZE;
            stop_pt = start_pt + (n * SIMSTEP_SIZE) - 1;

            if stop_pt > length(this.time)
                error("SIMULATION RAN OUT OF BOUNDS");
            end

            this.curr_step = this.curr_step + n;
            cut_time = this.time(start_pt:stop_pt);
            cut_carrier = this.carrier(start_pt:stop_pt);
            cut_data = this.data(start_pt:stop_pt);
            res = this.modulate_by_fsk(cut_time, cut_carrier, cut_data, this.channel);
        end

    end

end
