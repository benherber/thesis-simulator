% ------------------------ %
% Passive Tag Simulator    %
% Name   : Simulator.m     %
% Author : Benjamin Herber %
% Date   : Fall 2021       %
% ------------------------ %

classdef Simulator < handle
    %SIMULATOR object of backscatter tags.
    %   Simulator of system of backscatter tags communcating with a 
    %   basestation.

    properties (GetAccess = public, SetAccess = private)
        tags
        curr_step
        total_steps
        time
        carrier
        channel
        params
    end

    methods
        function this = Simulator(tags, tag_modes, channel, num_syms, sim_params)
            %SIMULATOR constructor
            %   Initialize and start simulator.

            this.params = sim_params;

            % Get time
            this.time = (0:(1 / this.params.Fs):(this.params.total_time - (1 / this.params.Fs)));
            
            % Carrier signal
            this.carrier = complex(this.params.amplitude * sin(complex(2 * pi * this.params.Fc * this.time)));

            this.channel = channel;

            % Init tags
            tag_sz = size(tags);
            this.tags = zeros(1, tag_sz(2));
            for idx = 1:numel(this.tags)
                % Get random data signal
                bits = randi([0, 1], num_syms, 1);
                data = repelem(bits, this.params.symb_sz).'; 

                this.tags(idx) = Tag(tags(1, idx), tags(2, idx), tags(3, idx), ...
                    tag_modes(idx), this.time, this.carrier, data, this.channel);
            end

            this.curr_step = 0;
            this.total_steps = num_syms * this.params.sim_symb_ratio;
        end

        function res = step(this)
            %STEP through one frame of the simulation
            %   Step one 1us frame in time through the simulation for the 
            %   entire system of tags.

            if this.curr_step > this.total_steps
                error("ATTEMPTED TO SIMULATE TOO MANY STEPS");
            end

            res = zeros(1, this.params.simstep_sz);

            for idx = 1:numel(this.tags)
                curr = this.tags(idx);
                res = res + curr.step();
            end

            this.curr_step = this.curr_step + 1;
        end
    end
end