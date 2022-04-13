clear; clc;

addpath(genpath(".."));
import bherber.thesis.*

consts_struct = struct( ...
    "Fc", 2.4e9, ...
    "Fb_base", 75e3, ...
    "Fb_step", 25e3, ...
    "Fb_channel_spacing", 50e3, ...
    "fs_granularity", 3, ...
    "num_symbs", 100, ...
    "symb_freq", 100e3, ...
    "sim_sym_ratio", 10, ...
    "amplitude", 0.1, ...
    "num_elements", 4 ...
    );
consts_cell = struct2cell(consts_struct);
simconsts = SimulationConstants(consts_cell{:});

% Get time
time = (0:(1 / simconsts.Fs):((simconsts.num_symbs * simconsts.symb_sz / simconsts.Fs) - (1 / simconsts.Fs)));

% Carrier signal
carrier = complex(simconsts.amplitude * sin(complex(2 * pi * simconsts.Fc * time)));

% Channel Definition
channel = @(given) awgn(given, 10, "measured");

powers = zeros(1, 100);

for jdx = 1:100
%     % Get random data signal
%     bits = randi([0, 1], 1, simconsts.num_symbs);
%     data = repelem(bits, simconsts.symb_sz); 
%     
%     % OOK Tag
%     
%     tag = Tag(0, 0, 0, TagType.FSK_LO, time, carrier, data, channel, simconsts);
%     
%     % Modulate
%     modded_steps = zeros(simconsts.simstep_sz, (simconsts.num_symbs * simconsts.sim_sym_ratio));
%     for idx = 1:(simconsts.num_symbs * simconsts.sim_sym_ratio)
%         modded_steps(:, idx) = tag.step();
%     end

    powers(jdx) = sum(abs(carrier(:)).^2)/numel(carrier); % rms(modded_steps(:)) ^ 2;
end

final_power = mean(powers);