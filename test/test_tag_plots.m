clear; clc;

addpath(genpath(".."));
import bherber.thesis.TagType

bits = [0,1,1,0,1,0,1,0,0,1,0];
data = repelem(bits, simconsts.symb_sz);
distance = 1;
num_steps = (simconsts.num_symbs + 1) * simconsts.sim_sym_ratio;
time = ((0:(double(simconsts.simstep_sz * (num_steps - simconsts.sim_sym_ratio)) - 1)) / simconsts.Fs);
delay_samples = @(distance) ceil((2 * distance / physconst("Lightspeed")) * simconsts.Fs);

%% OOK

ook_curr = bherber.thesis.tags.OOKTag(distance, 0, 0, bherber.thesis.TagType.OOK, bits, simconsts);
ook_modded_steps = zeros(simconsts.simstep_sz, num_steps);
for jdx = 1:num_steps
    ook_modded_steps(:, jdx) = ook_curr.step();
end
ook_mod = delayseq(ook_modded_steps(:), -1 * delay_samples(ook_curr.distance));
ook_mod = ook_mod(1:((num_steps - simconsts.sim_sym_ratio) * simconsts.simstep_sz)).';

figure;
hold on;

plot(time, ook_mod)
set(gca, "FontSize", 12);
xlabel("time(s)", FontSize=18);
ylabel("amplitude", FontSize=18);

xlines = zeros(1, simconsts.num_symbs);
for idx = 1:length(xlines)
    xlines(idx) = idx * (1 / simconsts.symb_freq);
end
xline(xlines, "--k");
xticks(xlines);
xlabels = strings(1, simconsts.num_symbs);
for idx = 1:length(xlabels)
    xlabels(idx) = string(idx) + "T_s";
end
xticklabels(xlabels);

%% MPAM

simconsts = bherber.thesis.SimulationConstants( ...
                2.4e9, ...   "Fc"
                1e6, ...   "Fb_base"
                1e6, ...   "Fb_step"
                0.5 / (1 / (100e3)), ...   "Fb_channel_spacing"
                1, ...       "num_channels"
                4, ...       "fs_granularity"
                10, ... "num_symbs"
                1e6, ...     "symb_freq"
                1, ...       "sim_sym_ratio"
                sqrt(8), ... "amplitude"
                4,  ...      "num_elements"
                3, ...        "pattern_len"
                4, ... "m_ary_modulation"
                1 ... "preamble_len"
            );

bits = [0,1,1,0,1,0,1,1,0,1,0,0,1,1,0,1,1,0,1,1,0,0];
ook_curr = bherber.thesis.tags.OOKTag(distance, 0, 0, bherber.thesis.TagType.OOK, bits, simconsts);
ook_modded_steps = zeros(simconsts.simstep_sz, num_steps);
for jdx = 1:num_steps
    ook_modded_steps(:, jdx) = ook_curr.step();
end
ook_mod = delayseq(ook_modded_steps(:), -1 * delay_samples(ook_curr.distance));
ook_mod = ook_mod(1:((num_steps - simconsts.sim_sym_ratio) * simconsts.simstep_sz)).';

figure;
hold on;

plot(time, ook_mod)
set(gca, "FontSize", 12);
xlabel("time(s)", FontSize=18);
ylabel("amplitude", FontSize=18);

xlines = zeros(1, simconsts.num_symbs);
for idx = 1:length(xlines)
    xlines(idx) = idx * (1 / simconsts.symb_freq);
end
xline(xlines, "--k");
xticks(xlines);
xlabels = strings(1, simconsts.num_symbs);
for idx = 1:length(xlabels)
    xlabels(idx) = string(idx) + "T_s";
end
xticklabels(xlabels);

%% FSK
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

bits = [0,1,1,0,1,0,1,0,0,1];
bits = [repelem(bits, 100), 0];
num_steps = (simconsts.num_symbs + 1) * simconsts.sim_sym_ratio;
% bits = randi([0,1], 1000, 1);
fsk_cur = bherber.thesis.tags.FSKTag(distance, 0, 0, bherber.thesis.TagType.FSK_LO, bits, simconsts, simconsts.freq_channels(1, :));
fsk_modded_steps = zeros(simconsts.simstep_sz, num_steps);
for jdx = 1:num_steps
    fsk_modded_steps(:, jdx) = fsk_cur.step();
end
fsk_mod = delayseq(fsk_modded_steps(:), -1 * delay_samples(fsk_cur.distance));
fsk_mod = fsk_mod(1:((num_steps - simconsts.sim_sym_ratio) * simconsts.simstep_sz)).';

nfft = 100000000;
f_fft = fftshift(fft(fsk_mod, nfft));
f = simconsts.Fs * (-nfft/2 : nfft/2 - 1) / nfft;
scaled_f = (f - 2.4e9) / 1e3;

f1 = figure();
hold on
scaled_fft = 2 * abs(f_fft) / length(fsk_mod);
plot(scaled_f, scaled_fft);
carrier_line = xline(0, "r", "Carrier", DisplayName="Carrier Frequency");
carrier_line.LabelVerticalAlignment = "middle";
carrier_line.LabelHorizontalAlignment = "center";

xlim([-2e3, 2e3])
legend("Modulated Signal", "Carrier Frequency","","");
ylabel("Amplitude", FontSize=18);
xlabel("Offset from Carrier Frequency (kHz)", FontSize=18);
% xticks(-2e3:50:500);

%% FSK - Spectrogram

simconsts = bherber.thesis.SimulationConstants( ...
                2.4e9, ...   "Fc"
                1e6, ...   "Fb_base"
                1e6, ...   "Fb_step"
                0.5 / (1 / (100e3)), ...   "Fb_channel_spacing"
                1, ...       "num_channels"
                4, ...       "fs_granularity"
                10, ... "num_symbs"
                1e6, ...     "symb_freq"
                1, ...       "sim_sym_ratio"
                sqrt(8), ... "amplitude"
                4,  ...      "num_elements"
                3, ...        "pattern_len"
                2, ... "m_ary_modulation"
                1 ... "preamble_len"
            );

bits = [0,1,1,0,1,0,1,0,0,1];
num_steps = (simconsts.num_symbs + 1) * simconsts.sim_sym_ratio;
fsk_cur = bherber.thesis.tags.FSKTag(distance, 0, 0, bherber.thesis.TagType.FSK_LO, bits, simconsts, simconsts.freq_channels(1, :));
fsk_modded_steps = zeros(simconsts.simstep_sz, num_steps);
for jdx = 1:num_steps
    fsk_modded_steps(:, jdx) = fsk_cur.step();
end
fsk_mod = delayseq(fsk_modded_steps(:), -1 * delay_samples(fsk_cur.distance));
fsk_mod = fsk_mod(1:((num_steps - simconsts.sim_sym_ratio) * simconsts.simstep_sz)).';

nfft = 100000;
f = figure();
spectrogram(fsk_mod, hamming(simconsts.symb_sz), ...
    0, nfft, simconsts.Fs, "yaxis", "power");
axs = f.Children;
data_objs = axs(2).Children;
% tmpYData = data_objs.YData;
% data_objs(1).YData = (tmpYData - 2.4);
ylim([2.4, 2.4025]);
yticks([0, 2.401, 2.4015,2.4025]);
yline([2.401, 2.4015])
yticklabels(string([0, 1000, 1500, 2500]));
clim([-80, -60])

xlines = zeros(1, simconsts.num_symbs);
for idx = 1:length(xlines)
    xlines(idx) = idx * (1 / simconsts.symb_freq) * 1e6;
end
xline(xlines, "--k");
xticks(xlines);
xlabels = strings(1, simconsts.num_symbs);
for idx = 1:length(xlabels)
    xlabels(idx) = string(idx) + "T_s";
end
xticklabels(xlabels);

ylabel("Offset from Carrier Frequency (kHz)");
xlabel("Symbol Time")

%% <FSK
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
                4, ... "m_ary_modulation"
                1 ... "preamble_len"
            );

bits = [0,0,1,1,0,1,1,0];
bits = repelem(bits.', 1, 1000/4);
bits = [bits(:).', 0, 0];
num_steps = (simconsts.num_symbs + 1) * simconsts.sim_sym_ratio;
% bits = randi([0,1], 1000, 1);
fsk_cur = bherber.thesis.tags.FSKTag(distance, 0, 0, bherber.thesis.TagType.FSK_LO, bits, simconsts, simconsts.freq_channels(1, :));
fsk_modded_steps = zeros(simconsts.simstep_sz, num_steps);
for jdx = 1:num_steps
    fsk_modded_steps(:, jdx) = fsk_cur.step();
end
fsk_mod = delayseq(fsk_modded_steps(:), -1 * delay_samples(fsk_cur.distance));
fsk_mod = fsk_mod(1:((num_steps - simconsts.sim_sym_ratio) * simconsts.simstep_sz)).';

nfft = 100000000;
f_fft = fftshift(fft(fsk_mod, nfft));
f = simconsts.Fs * (-nfft/2 : nfft/2 - 1) / nfft;
scaled_f = (f - 2.4e9) / 1e3;

f1 = figure();
hold on
scaled_fft = 2 * abs(f_fft) / length(fsk_mod);
plot(scaled_f, scaled_fft);
carrier_line = xline(0, "r", "Carrier", DisplayName="Carrier Frequency");
carrier_line.LabelVerticalAlignment = "middle";
carrier_line.LabelHorizontalAlignment = "center";

xlim([-3e3, 3e3])
legend("Modulated Signal", "Carrier Frequency","","");
ylabel("Amplitude", FontSize=18);
xlabel("Offset from Carrier Frequency (kHz)", FontSize=18);
% xticks(-2e3:50:500);

%% MFSK - Spectrogram

simconsts = bherber.thesis.SimulationConstants( ...
                2.4e9, ...   "Fc"
                1e6, ...   "Fb_base"
                1e6, ...   "Fb_step"
                0.5 / (1 / (100e3)), ...   "Fb_channel_spacing"
                1, ...       "num_channels"
                4, ...       "fs_granularity"
                10, ... "num_symbs"
                1e6, ...     "symb_freq"
                1, ...       "sim_sym_ratio"
                sqrt(8), ... "amplitude"
                4,  ...      "num_elements"
                3, ...        "pattern_len"
                4, ... "m_ary_modulation"
                1 ... "preamble_len"
            );

bits = [0,0,1,1,0,1,1,0,1,1,1,0,0,1,0,0,1,0,1,0,0,0];
num_steps = (simconsts.num_symbs + 1) * simconsts.sim_sym_ratio;
fsk_cur = bherber.thesis.tags.FSKTag(distance, 0, 0, bherber.thesis.TagType.FSK_LO, bits, simconsts, simconsts.freq_channels(1, :));
fsk_modded_steps = zeros(simconsts.simstep_sz, num_steps);
for jdx = 1:num_steps
    fsk_modded_steps(:, jdx) = fsk_cur.step();
end
fsk_mod = delayseq(fsk_modded_steps(:), -1 * delay_samples(fsk_cur.distance));
fsk_mod = fsk_mod(1:((num_steps - simconsts.sim_sym_ratio) * simconsts.simstep_sz)).';

nfft = 100000;
f = figure();
spectrogram(fsk_mod, hamming(simconsts.symb_sz), ...
    0, nfft, simconsts.Fs, "yaxis", "power");
axs = f.Children;
data_objs = axs(2).Children;
% tmpYData = data_objs.YData;
% data_objs(1).YData = (tmpYData - 2.4);
ylim([2.4, 2.4025]);
yticks([0, 2.401, 2.4011667, 2.401333, 2.4015, 2.4025]);
yline([2.401, 2.4011667, 2.401333, 2.4015])
yticklabels(string([0, 1000, 1167, 1333, 1500, 2500]));
clim([-80, -60])

xlines = zeros(1, simconsts.num_symbs);
for idx = 1:length(xlines)
    xlines(idx) = idx * (1 / simconsts.symb_freq) * 1e6;
end
xline(xlines, "--k");
xticks(xlines);
xlabels = strings(1, simconsts.num_symbs);
for idx = 1:length(xlabels)
    xlabels(idx) = string(idx) + "T_s";
end
xticklabels(xlabels);

ylabel("Offset from Carrier Frequency (kHz)");
xlabel("Symbol Time")

