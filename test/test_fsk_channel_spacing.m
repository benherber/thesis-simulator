clear; clc;

addpath(genpath(".."));
import bherber.thesis.*

% Possible spacings
spacings = (0:1e3:100e3);
CHANNEL0 = 0;
CHANNEL1 = 1;
bers = zeros(2, length(spacings));

parfor idx = 1:length(spacings)

    % Simulation Constants
    simconsts = params(spacings(idx));

    % Get time
    time = (0:(1 / simconsts.Fs):(simconsts.total_time - (1 / simconsts.Fs)));
    
    % Carrier signal
    carrier = complex(simconsts.amplitude * sin(complex(2 * pi * simconsts.Fc * time)));
    
    % Channel Definition
    channel = @(given) given;
    
    % Get random data signals
    bits1 = randi([0, 1], 1, simconsts.num_symbs);
    bits2 = randi([0, 1], 1, simconsts.num_symbs);
    data1 = repelem(bits1, simconsts.symb_sz);
    data2 = repelem(bits2, simconsts.symb_sz);

    % Simulate Tag Response
    tag1 = Tag(0, 0, 0, TagType.FSK_LO, time, carrier, data1, channel, simconsts);
    tag2 = Tag(0, 0, 0, TagType.FSK_HI, time, carrier, data2, channel, simconsts);

    % Modulate
    modded_steps = zeros(simconsts.simstep_sz, (simconsts.num_symbs * simconsts.sim_sym_ratio));
    for jdx = 1:(simconsts.num_symbs * simconsts.sim_sym_ratio)
        modded_steps(:, jdx) = tag2.step() + tag1.step();
    end
    f_modulation = modded_steps(:);
    
    % Demodulate
    modded_symbs = reshape(modded_steps, [simconsts.symb_sz, simconsts.num_symbs]);
    
    res_bits1 = zeros(1, numel(bits1));
    res_bits2 = zeros(1, numel(bits2));
    symb_times = reshape(time, [simconsts.symb_sz, simconsts.num_symbs]);
    carrier_split = reshape(carrier, [simconsts.symb_sz, simconsts.num_symbs]);
    channel0f1 = simconsts.fsk_channel0.f1;
    channel0f0 = simconsts.fsk_channel0.f0;
    channel1f1 = simconsts.fsk_channel1.f1;
    channel1f0 = simconsts.fsk_channel1.f0;
    for jdx = 1:simconsts.num_symbs
        res_bits1(jdx) = Tag.fsk_demodulate(modded_symbs(:, jdx), carrier_split(:, jdx), ... 
            symb_times(:, jdx), channel0f1, channel0f0);
        res_bits2(jdx) = Tag.fsk_demodulate(modded_symbs(:, jdx), carrier_split(:, jdx), ... 
            symb_times(:, jdx), channel1f1, channel1f0);
    end

    bers(:, idx) = [100 * (biterr(bits1, res_bits1) / simconsts.num_symbs), ...
        100 * (biterr(bits2, res_bits2) / simconsts.num_symbs)];

end

%% Plot Dependence

fig = figure();
hold on
plot(spacings / 1e3, bers)
title("Bit-Error Response Due to Inter-Channel Spacing")
ylabel("BER (%)");
xlabel("Inter-Channel Spacing (kHz)")
legend(["Low Frequency", "High Frequency"])
savefig(sprintf("plots/channel_spacing_fig::%s.fig", datestr(now, "mm-dd-yy-HH:MM:SS")));
hold off
close(fig);

%% Local Functions

function simconsts = params(spacing)
    import bherber.thesis.SimulationConstants
    consts_struct = struct( ...
        "Fc", 2.4e9, ...
        "Fb_base", 75e3, ...
        "Fb_step", 25e3, ...
        "Fb_channel_spacing", spacing, ...
        "fs_granularity", 3, ...
        "num_symbs", 1000, ...
        "symb_freq", 100e3, ...
        "sim_sym_ratio", 10, ...
        "amplitude", 0.1, ...
        "num_elements", 4 ...
        );
    consts_cell = struct2cell(consts_struct);
    simconsts = SimulationConstants(consts_cell{:});
end