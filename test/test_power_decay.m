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

power = @(signal) rms(signal)^2;

%% Power Calculations

MAX_DISTANCE = 100;
distances = 1:1:MAX_DISTANCE;
powers = zeros(size(distances));

num_steps = simconsts.num_symbs * simconsts.sim_sym_ratio;
step_sz = simconsts.simstep_sz;
parfor idx = 1:length(distances)
    curr = Tag(distances(idx), 0, 0, TagType.FSK_LO, time, carrier, data, channel, simconsts);
    modded_steps = zeros(step_sz, num_steps);
    for jdx = 1:num_steps
        modded_steps(:, jdx) = curr.step();
    end  
    powers(idx) = power(modded_steps(:));
end

%% Plot

fig = figure("Renderer", "painters", "Position", [10 10 700 600]);
semilogy(distances, powers);
title("Received Backscattered Power Relation to Tag Distance", FontSize=20);
ylabel("Received Power (rms)", FontSize=18);
xlabel("Tag Distance (m)", FontSize=18);
grid on
savefig(sprintf("plots/power_decay_fig::%s.fig", datestr(now, "mm-dd-yy-HH:MM:SS")));
close(fig);

            