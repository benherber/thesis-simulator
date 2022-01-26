classdef SimulationConstants < handle
    %SIMULATION_CONSTANTS defining model of system.
    %   All constants needed to run a simulation of a basestation with
    %   multiple backscatter tags in a system.

    properties (GetAccess = public, SetAccess = private)
        Fc % Carrier Frequency
        fsk_channel0 % First frequency channel for fsk
        fsk_channel1 % Second frequency channel for fsk
        wavelen {mustBePositive} % Wavelength of Carrier
        fs_granularity {mustBePositive} % Sampling Granularity
        Fs {mustBePositive} % Sampling frequency
        num_symbs {mustBePositive, mustBeGreaterThanOrEqual(num_symbs, 8)} % Total number of symbols in simulation
        symb_freq {mustBePositive} % Frequency of symbols
        symb_sz {mustBePositive} % Samples per symbol
        sim_sym_ratio {mustBePositive} % Simulation steps per symbol window
        simstep_freq {mustBePositive} % Frequency of simulation-steps
        simstep_sz  {mustBePositive}% Number of samples per simulation-step
        total_time {mustBePositive} % Total simulation time
        total_samples {mustBePositive} % Total simulation samples
        amplitude {mustBePositive} % Amplitude of Carrier
        num_elements {mustBePositive} % Number of elements in Van Atta Array
    end

    methods

        function this = SimulationConstants( ...
                Fc, ...
                Fb_base, ...
                Fb_step, ...
                Fb_channel_spacing, ...
                fs_granularity, ...
                num_symbs, ...
                symb_freq, ...
                sim_sym_ratio, ...
                amplitude, ...
                num_elements ...
            )
            %SIMULATION_CONSTANTS Construct an instance of constants
            %   Define constants for the simulation

            arguments
                Fc
                Fb_base
                Fb_step
                Fb_channel_spacing
                fs_granularity
                num_symbs
                symb_freq
                sim_sym_ratio
                amplitude
                num_elements
            end

            this.Fc = Fc;
            this.fsk_channel0 = struct("f0", Fb_base, "f1", Fb_base + Fb_step);
            this.fsk_channel1 = struct("f0", Fb_base + Fb_step + Fb_channel_spacing, ...
                "f1", Fb_base + Fb_step + Fb_channel_spacing + Fb_step);
            this.wavelen = physconst("Lightspeed") / Fc;
            this.fs_granularity = fs_granularity;
            this.Fs = Fc * fs_granularity;
            this.num_symbs = num_symbs;
            this.symb_freq = symb_freq;
            this.symb_sz = (1 / symb_freq) * this.Fs;
            this.sim_sym_ratio = sim_sym_ratio;
            this.simstep_freq = symb_freq * sim_sym_ratio;
            this.simstep_sz = (1 / (symb_freq * sim_sym_ratio)) * this.Fs;
            this.total_time = num_symbs * (1 / symb_freq);
            this.total_samples = this.symb_sz * num_symbs;
            this.amplitude = amplitude;
            this.num_elements = num_elements;

        end

    end

end
