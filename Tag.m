% ------------------------ %
% Back Scatter Tag         %
% Name   : Tag.m           %
% Author : Benjamin Herber %
% Date   : Fall 2021       %
% ------------------------ %

classdef Tag < handle
    %TAG model
    %   Backscatter tag model of a passive Van Atta configuration

%%  Public Properties

    properties
        x % x-position of tag relative to basestation
        y % y-position of tag relative to basestation
        z % z-position of tag relative to basestation
        mode % Frequency mode ('lo', 'hi', or 'random' bands)
    end

%%  Private Properties

    properties (GetAccess = public, SetAccess = private)
        curr_step % Current time-slice in modulation
        carrier % Delayed carrier signal
        data % Delayed data signal
        time % Time signal
        channel % Propagation channel
        params % Simulation parameters
    end

%%  Dependent Properties

    properties (Dependent)
        distance % Tag distance from basestation
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

            switch this.mode
                case "lo"
                    f1 = this.params.fsk_channel0.f1;
                    f0 = this.params.fsk_channel0.f0;

                case "hi"
                    f1 = this.params.fsk_channel1.f1;
                    f0 = this.params.fsk_channel1.f0;

                case "random"

                    if logical(randi([0, 1]))
                        f1 = this.params.fsk_channel0.f1;
                        f0 = this.params.fsk_channel0.f0;
                    else
                        f1 = this.params.fsk_channel1.f1;
                        f0 = this.params.fsk_channel1.f0;
                    end

                otherwise
                    error("UNSUPPORTED MODE");
            end

        end

% ----------------------------------------------------------------------- %

        function modulated = modulate_by_fsk(this, t, carrier, data, channel)
            % MODULATE_BY_FSK a given signal.
            %   Given carrier and data signals, modulate the carrier wave according to
            %   frequency-shift keying as a backscatter tag would through a given channel.

            % Frequency modulation

            % 1. Add channel noise headed to tag
            noisy_carrier = carrier;%friis_path_loss(channel(carrier), this.distance, this.params);

            % 2. Modulate
            f_modulation = zeros(size(data));
            [f1, f0] = this.freq();
            sq_one = square((f1) * 2 * pi * t);
            sq_zero = square((f0) * 2 * pi * t);
            all_ones = sq_one .* noisy_carrier * ...
                vanatta_gain(this.params.num_elements, this.params.Fc, f1);
            all_zeros = sq_zero .* noisy_carrier * ...
                vanatta_gain(this.params.num_elements, this.params.Fc, f0);

            for idx = 1:length(f_modulation)

                if data(idx) == 1
                    f_modulation(idx) = all_ones(idx);
                else
                    f_modulation(idx) = all_zeros(idx);
                end

            end

            % 3. Add channel noise headed back to basestation
            modulated = f_modulation; %friis_path_loss(channel(f_modulation), this.distance, this.params);
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
            this.curr_step = 0;
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

            if nargin < 2
                n = 1;
            end

            start_pt = uint64(this.curr_step * this.params.simstep_sz + 1);
            stop_pt = uint64(start_pt + (n * this.params.simstep_sz) - 1);

            if stop_pt > length(this.time)
                error("SIMULATION RAN OUT OF BOUNDS");
            end

            cut_time = this.time(start_pt:stop_pt);
            cut_carrier = this.carrier(start_pt:stop_pt);
            cut_data = this.data(start_pt:stop_pt);
            res = this.modulate_by_fsk(cut_time, cut_carrier, cut_data, this.channel);
            this.curr_step = this.curr_step + n;
        end

    end

end
