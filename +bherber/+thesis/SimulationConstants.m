classdef SimulationConstants < handle
    %SIMULATION_CONSTANTS defining model of system.
    %   All constants needed to run a simulation of a basestation with
    %   multiple backscatter tags in a system.

    properties (GetAccess = public, SetAccess = immutable)
        Fc % Carrier Frequency
        freq_channels % Frequency channels for fsk/FSK-FH
        wavelen {mustBePositive} % Wavelength of Carrier
        fs_granularity {mustBePositive} % Sampling Granularity
        interchannel_spacing {mustBePositive} % Spacing between frequency channels
        pattern_len {mustBePositive} % Length of Freq Hop Pattern
        num_channels {mustBePositive} % Number of frequency channels available
        Fs {mustBePositive} % Sampling frequency
        num_symbs {mustBePositive} % Total number of symbols in simulation
        symb_freq {mustBePositive} % Frequency of symbols
        symb_sz {mustBePositive} % Samples per symbol
        sim_sym_ratio {mustBePositive} % Simulation steps per symbol window
        simstep_freq {mustBePositive} % Frequency of simulation-steps
        simstep_sz  {mustBePositive}% Number of samples per simulation-step
        total_time {mustBePositive} % Total simulation time
        total_samples {mustBePositive} % Total simulation samples
        amplitude {mustBePositive} % Amplitude of Carrier
        num_elements {mustBePositive} % Number of elements in Van Atta Array
        m_ary_modulation {mustBePositive, mustBeGreaterThanOrEqual(m_ary_modulation, 2)} % M-Ary Mod. Scheme
        m_ary_amplitudes % Amplitudes for M-Ary ASK
        slots_per_frame {mustBePositive} % Number of TDM slots in a given frame
        preamble_len {mustBePositive} % Length in bits of preamble ids
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
                pattern_len, ...
                m_ary_scheme, ...
                preamble_len ...
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
                m_ary_scheme = 2
                preamble_len = 16
            end

            this.Fc = Fc;
%             this.freq_channels = [struct("f0", 4e6, "f1", 8e6)];
            this.freq_channels = zeros(num_channels, m_ary_scheme);
            freq_spacing = 1 / (2 * (1 / symb_freq));
            bands = linspace(0, freq_spacing, m_ary_scheme);
            for idx = 1:num_channels
                this.freq_channels(idx, :) = Fb_base + bands + ...
                    ((idx - 1) * (Fb_channel_spacing + bands(end)));
%                 f0 = Fb_base + ((idx - 1) * (Fb_step + Fb_channel_spacing));
%                 f1 = Fb_base + (idx * Fb_step) + ((idx - 1) * Fb_channel_spacing);
%                 this.freq_channels = [this.freq_channels, struct("f0", f0, "f1", f1)];
            end
            this.interchannel_spacing = Fb_channel_spacing;
            this.num_channels = num_channels;
            this.pattern_len = pattern_len;
            this.wavelen = physconst("Lightspeed") / Fc;
            this.fs_granularity = fs_granularity;
            this.Fs = Fc * fs_granularity;
            this.num_symbs = num_symbs;
            this.symb_freq = symb_freq;
            this.symb_sz = int32((1 / symb_freq) * this.Fs);
            this.sim_sym_ratio = sim_sym_ratio;
            this.simstep_freq = symb_freq * sim_sym_ratio;
            this.simstep_sz = (1 / (symb_freq * sim_sym_ratio)) * this.Fs;
            this.total_time = (num_symbs + 1) * (1 / symb_freq);
            this.total_samples = this.symb_sz * num_symbs;
            this.amplitude = amplitude;
            this.num_elements = num_elements;
            this.m_ary_modulation = m_ary_scheme;
            this.m_ary_amplitudes = linspace(0, amplitude, m_ary_scheme);
            this.preamble_len = preamble_len;
        end

    end

end
