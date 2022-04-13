clear; clc

import bherber.thesis.PacketSimulator ...
    bherber.thesis.SimulationConstants ...
    bherber.thesis.tags.FreqHopTag ...
    bherber.thesis.TagType


params = SimulationConstants( ...
                2.4e9, ...   "Fc"
                1e6, ...   "Fb_base"
                1e6, ...   "Fb_step"
                0.5 / (1 / (1e6)), ...   "Fb_channel_spacing"
                5, ...       "num_channels"
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

bits = repelem(0, 11);
hopset = [5,1,2,5,1,2,5,1,2,5,1,2];
tag = bherber.thesis.tags.FreqHopTag(1, 0, 0, TagType.FREQ_HOP, bits, params, params.freq_channels, hopset);
delay = 2 * tag.distance / physconst("Lightspeed");
delay_samples = ceil(delay * params.Fs);

steps = zeros(params.symb_sz, params.num_symbs + 1);
for idx = 1:(params.num_symbs + 1)
    steps(:, idx) = tag.step();
end

f_modulation = steps(:);
f_modulation = delayseq(f_modulation, -1 * delay_samples);
f_modulation = f_modulation(1:(params.num_symbs * params.symb_sz));

% %%
% 
% symbs = f_modulation;
% 
% nfft = 100000000;
% f_fft = fftshift(fft(f_modulation, nfft));
% 
% % Frequency
% f = simconsts.Fs * (-nfft/2 : nfft/2 - 1) / nfft;
% scaled_f = (f - 2.4e9) / 1e3;
% 
% %% FFT Plot
% 
% f1 = figure();
% hold on
% scaled_fft = 2 * abs(f_fft) / length(f_modulation);
% plot(scaled_f, scaled_fft);

%% 

nfft = 100000;
f = figure();
spectrogram(f_modulation, hamming(params.symb_sz), ...
    0, nfft, params.Fs, "yaxis", "power");
axs = f.Children;
data_objs = axs(2).Children;
tmpYData = data_objs.YData;
% data_objs(1).YData = (tmpYData - 2.4);
ylim([2.4, 2.406]);
clim([-80, -60])

yticks([2.4, 2.401, 2.402, 2.403, 2.404, 2.405, 2.406]);
yline([2.401, 2.402, 2.403, 2.404, 2.405])
yticklabels(string([0, 1000, 2000, 3000, 4000, 5000, 6000]));

xlines = zeros(1, params.num_symbs);
for idx = 1:length(xlines)
    xlines(idx) = idx * (1 / params.symb_freq) * 1e6;
end
xline(xlines, "--k");
xticks(xlines);
xlabels = strings(1, params.num_symbs);
for idx = 1:length(xlabels)
    xlabels(idx) = string(idx) + "T_s";
end
xticklabels(xlabels);

ylabel("Offset from Carrier Frequency (kHz)", FontSize=18);
xlabel("Symbol Time", FontSize=18)