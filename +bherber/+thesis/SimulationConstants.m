classdef SimulationConstants < handle
    %SIMULATION_CONSTANTS defining model of system.
    %   All constants needed to run a simulation of a basestation with
    %   multiple backscatter tags in a system.

    properties (GetAccess = public, SetAccess = private)
        Fc % Carrier Frequency
        freq_channels % Frequency channels for fsk/FSK-FH
        wavelen {mustBePositive} % Wavelength of Carrier
        fs_granularity {mustBePositive} % Sampling Granularity
        interchannel_spacing {mustBePositive} % Spacing between frequency channels
        pattern_len {mustBePositive} % Length of Freq Hop Pattern
        num_channels {mustBePositive} % Number of frequency channels available
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
                num_channels, ...
                fs_granularity, ...
                num_symbs, ...
                symb_freq, ...
                sim_sym_ratio, ...
                amplitude, ...
                num_elements, ...
                pattern_len ...
            )
            %SIMULATION_CONSTANTS Construct an instance of constants
            %   Define constants for the simulation

            arguments
                Fc
                Fb_base
                Fb_step
                Fb_channel_spacing
                num_channels
                fs_granularity
                num_symbs
                symb_freq
                sim_sym_ratio
                amplitude
                num_elements
                pattern_len = 0
            end

            this.Fc = Fc;
            this.freq_channels = [];
            for idx = 1:num_channels
                f0 = Fb_base + ((idx - 1) * (Fb_step + Fb_channel_spacing));
                f1 = Fb_base + (idx * Fb_step) + ((idx - 1) * Fb_channel_spacing);
                this.freq_channels = [this.freq_channels, struct("f0", f0, "f1", f1)];
            end
            this.interchannel_spacing = Fb_channel_spacing;
            this.num_channels = num_channels;
            this.pattern_len = pattern_len;
            this.wavelen = physconst("Lightspeed") / Fc;
            this.fs_granularity = fs_granularity;
            this.Fs = Fc * fs_granularity;
            this.num_symbs = num_symbs;
            this.symb_freq = symb_freq;
            this.symb_sz = (1 / symb_freq) * this.Fs;
            this.sim_sym_ratio = sim_sym_ratio;
            this.simstep_freq = symb_freq * sim_sym_ratio;
            this.simstep_sz = (1 / (symb_freq * sim_sym_ratio)) * this.Fs;
            this.total_time = (num_symbs + 1) * (1 / symb_freq);
            this.total_samples = this.symb_sz * num_symbs;
            this.amplitude = amplitude;
            this.num_elements = num_elements;

        end

    end

end
