clear; clc;

addpath(genpath(".."));
import bherber.thesis.*

consts_struct = struct( ...
    "Fc", 2.4e9, ...
    "Fb_base", 75e3, ...
    "Fb_step", 25e3, ...
    "Fb_channel_spacing", 25e3, ...
    "fs_granularity", 3, ...
    "num_symbs", 1000, ...
    "symb_freq", 100e3, ...
    "sim_sym_ratio", 10, ...
    "amplitude", 0.1, ...
    "num_elements", 4 ...
    );
consts_cell = struct2cell(consts_struct);
simconsts = SimulationConstants(consts_cell{:});

% Get time
time = (0:(1 / simconsts.Fs):(simconsts.total_time - (1 / simconsts.Fs)));

% Carrier signal
carrier = complex(simconsts.amplitude * sin(complex(2 * pi * simconsts.Fc * time)));

% Channel Definition
channel = @(given) given;

% Get random data signal
bits = randi([0, 1], simconsts.num_symbs, 1);
data = repelem(bits, simconsts.symb_sz);

%% Power Calculations

MAX_DISTANCE = 100;
distances = 1:1:MAX_DISTANCE;
bers = zeros(length(distances), 2);

num_steps = simconsts.num_symbs * simconsts.sim_sym_ratio;
step_sz = simconsts.simstep_sz;
parfor idx = 1:length(distances)
    fsk_curr = Tag(distances(idx), 0, 0, TagType.FSK_LO, time, carrier, data, channel, simconsts);
    fsk_modded_steps = zeros(step_sz, num_steps);

    ook_curr = Tag(distances(idx), 0, 0, TagType.OOK, time, carrier, data, channel, simconsts);
    ook_modded_steps = zeros(step_sz, num_steps);
    for jdx = 1:num_steps
        fsk_modded_steps(:, jdx) = fsk_curr.step();
        ook_modded_steps(:, jdx) = ook_curr.step();
    end  
    
    bers(idx, :) = [ook_ber(bits, time, ook_modded_steps, carrier, simconsts), ...
        fsk_ber(bits, time, fsk_modded_steps, carrier, TagType.FSK_LO, simconsts)];
end

%% Plot

fig = figure("Renderer", "painters", "Position", [10 10 700 600]);
plot(distances, bers)
title("Bit-Error Rate (BER) as a Function of Distance", FontSize=20);
ylabel("BER (%)", FontSize=18);
xlabel("Tag Distance (m)", FontSize=18);
legend(["OOK", "FSK"]);
grid on
savefig(sprintf("plots/distance_delay::%s.fig", datestr(now, "mm-dd-yy-HH:MM:SS")));
close(fig);

%% Calc BER

function res = fsk_ber(bits, time, modded_steps, carrier, tag_type, params)
    import bherber.thesis.TagType bherber.thesis.Tag

    modded_symbs = reshape(modded_steps, [params.symb_sz, params.num_symbs]);
    res_bits = zeros(numel(bits), 1);
    symb_times = reshape(time, [params.symb_sz, params.num_symbs]);
    carrier_split = reshape(carrier, [params.symb_sz, params.num_symbs]);

    if tag_type == TagType.FSK_HI
        f1 = params.fsk_channel1.f1;
        f0 = params.fsk_channel1.f0;
    else
        f1 = params.fsk_channel0.f1;
        f0 = params.fsk_channel0.f0;
    end

    for idx = 1:params.num_symbs
        res_bits(idx) = Tag.fsk_demodulate(modded_symbs(:, idx), carrier_split(:, idx), ... 
            symb_times(:, idx), f1, f0);
    end

    % Calculate BER
    res = 100 * (biterr(bits, res_bits) / params.num_symbs);
end

function res = ook_ber(bits, time, modded_steps, carrier, params)
    import bherber.thesis.TagType bherber.thesis.Tag

   modded_symbs = reshape(modded_steps, [params.symb_sz, params.num_symbs]);

    res_bits = zeros(numel(bits), 1);
    symb_times = reshape(time, [params.symb_sz, params.num_symbs]);
    carrier_split = reshape(carrier, [params.symb_sz, params.num_symbs]);
    for idx = 1:params.num_symbs
        res_bits(idx) = Tag.ook_demodulate(modded_symbs(:, idx), carrier_split(:, idx), ... 
             symb_times(:, idx));
    end

    % Calculate BER
    res = 100 * (biterr(bits, res_bits) / params.num_symbs);
end

            