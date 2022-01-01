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

    properties (Access = private)
        tags
        curr_step
        total_steps
        time
        carrier
        channel
    end

    methods
        function this = Simulator(tags, tag_modes, channel, num_syms)
            %SIMULATOR constructor
            %   Initialize and start simulator.

            simulation_constants

            % Get time
            this.time = (0:1/Fs:(TOTAL_TIME - 1/Fs));
            
            % Carrier signal
            this.carrier = complex(A * sin(complex(2 * pi * FC * this.time)));

            this.channel = channel;

            % Init tags
            tag_sz = size(tags);
            this.tags = zeros(1, tag_sz(2));
            for idx = 1:numel(this.tags)
                % Get random data signal
                bits = randi([0, 1], num_syms, 1);
                data = repelem(bits, SYMB_SIZE).'; 

                this.tags(idx) = Tag(tags(1, idx), tags(2, idx), tags(3, idx), ...
                    tag_modes(idx), this.time, this.carrier, data, this.channel);
            end

            this.curr_step = 0;
            this.total_steps = num_syms * SIM_SYMB_RATIO;
        end

        function res = step(this)
            %STEP through one frame of the simulation
            %   Step one 1us frame in time through the simulation for the 
            %   entire system of tags.

            if this.curr_step > this.total_steps
                error("ATTEMPTED TO SIMULATE TOO MANY STEPS");
            end

            simulation_constants
            res = zeros(1, Fs * SIMSTEP);

            for idx = 1:numel(this.tags)
                curr = this.tags(idx);
                res = res + curr.step();
            end

            this.curr_step = this.curr_step + 1;
        end
    end
end