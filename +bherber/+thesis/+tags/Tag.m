% ------------------------ %
% Back Scatter Tag         %
% Name   : Tag.m           %
% Author : Benjamin Herber %
% Date   : Spring 2021     %
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
        bits (:, 1) {mustBeNonmissing, mustBeInRange(bits, 0, 1)} % bitstream
        is_active logical = true; % Tag is currently active in simulation window
        params bherber.thesis.SimulationConstants % Simulation parameters
        snr_db double % SNR of noise if applicable;
        complex_noise logical % Apply complex noise
    end

    %%  Dependent Properties

    properties (Dependent)
        distance % Tag distance from basestation
    end

    %%  Static Methods
    methods (Static)
        function res = vanatta_gain(num_elements, fc, fb)
            %VANATTA_GAIN of a backscatter tag in a Van Atta Configuration.
            %   Calculate gain of a backscatter tag in a Van Atta
            %   Configuration for given input and output frequencies.

            tmp = 0;
            C = physconst("Lightspeed");
            lambda_c = C / fc;
            lambda_b = C / fb;
            spacing = lambda_c / 2;

            find_phase = @(wavelen, spacing, n) ((2.0 * pi) / wavelen) * spacing * n;

            for idx = 1:num_elements
                tmp = tmp + exp(1i * ...
                    (find_phase(lambda_b, spacing, (idx - 1)) ...
                    - find_phase(lambda_c, spacing, (idx - 1))));
            end

            res = tmp;
        end

        % ----------------------------------------------------------------------- %

        function res = friis_path_loss(signal, distance, params)
            %FRIIS_PATH_LOSS of signal through free space
            %   Using Friis Formula for path loss through free space, calculate
            %   the remaining wave.

            if distance == 0
                res = signal;
            else
                % Calculate power loss
                amplitude_loss = params.wavelen / (4 * pi * distance);

                % Calculate the resultant wave
                res = signal * amplitude_loss;
            end
        end

    end

    %% Public Instance Methods

    methods

        function this = Tag(x, y, z, mode, bits, sim_params, options)
            %TAG constructor
            %   Create a tag a specific distance from a basestation in 3D
            %   space.
            
            arguments
                x double {mustBeScalarOrEmpty} % x-position of tag relative to basestation
                y double {mustBeScalarOrEmpty} % y-position of tag relative to basestation
                z double {mustBeScalarOrEmpty} % z-position of tag relative to basestation
                mode bherber.thesis.TagType % Modulation mode
                bits (1, :) {mustBeNonmissing, mustBeInRange(bits, 0, 1)} % bitstream
                sim_params bherber.thesis.SimulationConstants % Simulation parameters
                options.snr_db double = NaN;
                options.complex_noise logical = false;
            end

            this.params = sim_params;
            this.x = x;
            this.y = y;
            this.z = z;
            this.mode = mode;
            this.curr_step = int64(0);
            this.bits = [bits, 0];
            this.snr_db = options.snr_db;
            this.complex_noise = options.complex_noise;

        end

% ----------------------------------------------------------------------- %

        function activate(this)
            this.is_active = true;
        end

% ----------------------------------------------------------------------- %

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
    end

% ----------------------------------------------------------------------- %

    methods (Access = protected)

        function res = step_data(this)

            if this.curr_step > (((this.params.num_symbs + 1) * this.params.sim_sym_ratio) - 1)
                error("SIMULATION RAN OUT OF BOUNDS");
            end

            delay = 2 * this.distance / physconst("Lightspeed");
            time_slice = (0:(1 / this.params.Fs):(this.params.simstep_sz * (1 / this.params.Fs) - (1 / this.params.Fs))) ...
                + double(this.curr_step * this.params.simstep_sz) * (1 / this.params.Fs) ...
                - delay;
            invalid_times = time_slice < 0;

            carrier_slice = this.params.amplitude * cos(2 * pi * this.params.Fc * time_slice);
            carrier_slice(invalid_times) = 0;
            if ~isnan(this.snr_db)
                linear_snr = 10 ^ (this.snr_db / 10);
                linear_carrier_power = (this.params.amplitude ^ 2) / 2;
                carrier_noise = bherber.thesis.channel_models.AWGNChannelModel.noise(...
                    length(time_slice), linear_carrier_power, linear_snr, this.complex_noise);
                carrier_slice = carrier_slice + carrier_noise;
            end

            curr_symb = ceil((double(this.curr_step) + 1) / this.params.sim_sym_ratio);
            prev_steps_symb = ceil(double(this.curr_step) / this.params.sim_sym_ratio);
            if this.curr_step == 0
                data = repelem([0, this.bits(1)], this.params.simstep_sz);
            else
                bits_o_interest = [this.bits(prev_steps_symb), this.bits(curr_symb)];
                data = repelem(bits_o_interest, this.params.simstep_sz);
            end
            delay_samples = ceil(delay * this.params.Fs);
            data_start = this.params.simstep_sz + 1 - delay_samples;
            data_stop = data_start + this.params.simstep_sz - 1;
            data = data(data_start:data_stop);


            res = struct(...
                    "delay_samples", delay_samples, ...
                    "time_slice", time_slice, ...
                    "carrier_slice", carrier_slice, ...
                    "curr_symb", curr_symb, ...
                    "prev_steps_symb", prev_steps_symb, ...
                    "data", data ...
                );

            this.curr_step = this.curr_step + 1;
        end
    end

% ----------------------------------------------------------------------- %

    methods (Abstract)
        step(this)
        demodulate(this, signal)
    end

end
