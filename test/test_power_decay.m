clear; clc;

addpath(genpath(".."));
import bherber.thesis.TagType

simconsts = bherber.thesis.SimulationConstants( ...
                2.4e9, ...   "Fc"
                1e6, ...   "Fb_base"
                1e6, ...   "Fb_step"
                0.5 / (1 / (100e3)), ...   "Fb_channel_spacing"
                1, ...       "num_channels"
                4, ...       "fs_granularity"
                1000, ... "num_symbs"
                1e6, ...     "symb_freq"
                1, ...       "sim_sym_ratio"
                sqrt(8), ... "amplitude"
                4,  ...      "num_elements"
                3, ...        "pattern_len"
                2, ... "m_ary_modulation"
                1 ... "preamble_len"
            );

% Get random data signal
bits = randi([0, 1], simconsts.num_symbs, 1);
data = repelem(bits, simconsts.symb_sz);

power = @(signal) rms(signal)^2;

%% Power Calculations

MAX_DISTANCE = 100;
distances = 1:1:MAX_DISTANCE;
powers = zeros(length(distances), 2);

num_steps = simconsts.num_symbs * simconsts.sim_sym_ratio;
step_sz = simconsts.simstep_sz;
parfor idx = 1:length(distances)
    fsk_curr = bherber.thesis.tags.FSKTag(distances(idx), 0, 0, TagType.FSK_LO, bits, simconsts, simconsts.freq_channels(1, :));
    fsk_modded_steps = zeros(step_sz, num_steps);
    ook_curr = bherber.thesis.tags.OOKTag(distances(idx), 0, 0, TagType.OOK, bits, simconsts);
    ook_modded_steps = zeros(step_sz, num_steps);
    for jdx = 1:num_steps
        fsk_modded_steps(:, jdx) = fsk_curr.step();
        ook_modded_steps(:, jdx) = ook_curr.step();
    end  
    powers(idx, :) = [power(ook_modded_steps(:)), power(fsk_modded_steps(:))];
end

%% Plot

fig = figure("Renderer", "painters", "Position", [10 10 700 600]);
plot(distances, 10*log10(powers / 1e-3), LineWidth=3);
set(gca, "FontSize", 12);
% title("Received Backscattered Power Relation to Tag Distance", FontSize=20);
ylabel("Received Power (dBm)", FontSize=18);
xlabel("Tag Distance (m)", FontSize=18);
legend(["OOK", "FSK"])
grid on
% savefig(sprintf("plots/power_decay_fig::%s.fig", datestr(now, "mm-dd-yy-HH:MM:SS")));
% close(fig);

            