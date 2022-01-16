classdef SimulationConstants < handle
    %SIMULATION_CONSTANTS defining model of system.
    %   All constants needed to run a simulation of a basestation with
    %   multiple backscatter tags in a system.

    properties (GetAccess = public, SetAccess = private)
        Fc % Carrier Frequency
        fsk_channel0; % First frequency channel for fsk
        fsk_channel1; % Second frequency channel for fsk
        wavelen % Wavelength of Carrier
        fs_granularity % Sampling Granularity
        Fs % Sampling frequency
        num_symbs % Total number of symbols in simulation
        symb_freq % Frequency of symbols
        symb_sz % Samples per symbol
        sim_sym_ratio % Simulation steps per symbol window
        simstep_freq % Frequency of simulation-steps
        simstep_sz % Number of samples per simulation-step
        total_time % Total simulation time
        total_samples % Total simulation samples
        amplitude % Amplitude of Carrier
        num_elements % Number of elements in Van Atta Array
    end

    methods

        function this = SimulationConstants( ...
                Fc, ...
                Fb_base, ...
                Fb_step, ...
                fs_granularity, ...
                num_symbs, ...
                symb_freq, ...
                sim_sym_ratio, ...
                amplitude, ...
                num_elements ...
            )
            %SIMULATION_CONSTANTS Construct an instance of constants
            %   Define constants for the simulation

            this.Fc = Fc;
            this.fsk_channel0 = struct("f0", Fb_base, "f1", Fb_base + 1 * Fb_step);
            this.fsk_channel1 = struct("f0", Fb_base + 4 * Fb_step, "f1", Fb_base + 5 * Fb_step);
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
