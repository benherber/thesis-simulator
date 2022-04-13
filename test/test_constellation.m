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
time = (0:(1 / simconsts.Fs):(simconsts.total_time - (1 / simconsts.Fs)));

% Carrier signal
carrier = complex(simconsts.amplitude * sin(complex(2 * pi * simconsts.Fc * time)));

% Channel Definition
channel = @(given) awgn(given, 2, "measured");

% Get random data signal
bits = randi([0, 1], 1, simconsts.num_symbs);
data = repelem(bits, simconsts.symb_sz);

%% Run steps for modulation -- OOK

fig1 = figure();
hold on

tag = Tag(2, 0, 0, TagType.OOK, time, carrier, data, channel, simconsts);

% Modulate
modded_steps = zeros(simconsts.simstep_sz, (simconsts.num_symbs * simconsts.sim_sym_ratio));
for idx = 1:(simconsts.num_symbs * simconsts.sim_sym_ratio)
    modded_steps(:, idx) = tag.step();
end

% Demodulate
modded_symbs = reshape(modded_steps, [simconsts.symb_sz, simconsts.num_symbs]);

res_bits2 = zeros(1, numel(bits));
symb_times = reshape(time, [simconsts.symb_sz, simconsts.num_symbs]);
carrier_split = reshape(carrier, [simconsts.symb_sz, simconsts.num_symbs]);
parfor idx = 1:simconsts.num_symbs
    res_bits2(idx) = Tag.ook_demodulate(modded_symbs(:, idx), carrier_split(:, idx), ... 
         symb_times(:, idx));
end

% title(sprintf("Constellation Plot of OOK-Modulated Signal with BER=%%%0.2f", ...
%     (biterr(bits, res_bits2) / simconsts.num_symbs) * 100));

title("10 OOK-Modulated Symbols with Noise")

ook_mod = modded_symbs(:);
plot(time((1):(simconsts.symb_sz * 10)), real(ook_mod((1):(simconsts.symb_sz * 10))));
plot(time((1):(simconsts.symb_sz * 10)), imag(ook_mod((1):(simconsts.symb_sz * 10))));

savefig("ook_modulated_plot_wnoise.fig");

close(fig1);

%% Run steps for modulation -- FSK


fig2 = figure();
hold on

tag = Tag(2, 0, 0, TagType.FSK_LO, time, carrier, data, channel, simconsts);

% Modulate
modded_steps = zeros(simconsts.simstep_sz, (simconsts.num_symbs * simconsts.sim_sym_ratio));
for idx = 1:(simconsts.num_symbs * simconsts.sim_sym_ratio)
    modded_steps(:, idx) = tag.step();
end

% Demodulate
modded_symbs = reshape(modded_steps, [simconsts.symb_sz, simconsts.num_symbs]);

res_bits2 = zeros(1, numel(bits));
symb_times = reshape(time, [simconsts.symb_sz, simconsts.num_symbs]);
carrier_split = reshape(carrier, [simconsts.symb_sz, simconsts.num_symbs]);
parfor idx = 1:simconsts.num_symbs
    res_bits2(idx) = Tag.ook_demodulate(modded_symbs(:, idx), carrier_split(:, idx), ... 
         symb_times(:, idx));
end

% title(sprintf("Constellation Plot of FSK-Modulated Signal with BER=%%%0.2f", ...
%     (biterr(bits, res_bits2) / simconsts.num_symbs) * 100));

% savefig("FSKConstellation.fig");

title("10 FSK-Modulated Symbols with Noise")

fsk_mod = modded_symbs(:);
plot(time((1):(simconsts.symb_sz * 10)), real(fsk_mod((1):(simconsts.symb_sz * 10))));
plot(time((1):(simconsts.symb_sz * 10)), imag(fsk_mod((1):(simconsts.symb_sz * 10))));

savefig("fsk_modulated_plot_wnoise.fig");

close(fig2);