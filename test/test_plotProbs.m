clear; clc;
addpath(genpath(".."));
folder_name = "ProbSlurm";

%% Aggregate FTDMA Collision

files = dir(fullfile("data/" + folder_name + "/", "freqtdmcollisions*.mat"));
for idx = 1:length(files)
    load(fullfile(files(idx).folder, files(idx).name));
end

num_channels = 340;


num_tags = ceil(linspace(1, num_channels, 50));
num_tags = num_tags(1:17);
prob = zeros(size(num_tags));
for idx = 1:length(prob)
    prob(idx) = eval(sprintf("tags%dnumchannels340numslots10", num_tags(idx)));
end

figure();
plot(num_tags/(num_channels * 10), prob, LineWidth=3)
title("Collisions Given Maximum Parameter Constraints", FontSize=20);
xlabel("Activity Factor (\alpha)", FontSize=18);
ylabel("Collision (%)", FontSize=18);

%% Aggregate FTDMA Collision

clear -except folder_name; clc;
folder_name = "ProbSlurm";

files = dir(fullfile("data/" + folder_name + "/", "freqcollisions*.mat"));
for idx = 1:length(files)
    load(fullfile(files(idx).folder, files(idx).name));
end

prob50 = zeros(1, 50);
tags50 = linspace(1, 50, 50);
for idx = 1:length(prob50)
    prob50(idx) = eval(sprintf("tags%dnumchannels50", tags50(idx)));
end

tags350 = ceil(linspace(1, 350, 50));
tags350 = tags350(tags350 <= 330);
prob350 = zeros(size(tags350));
for idx = 1:length(tags350)
    exp = sprintf("tags%dnumchannels350", tags350(idx));
    prob350(idx) = eval(exp);
end

tags500 = ceil(linspace(1, 500, 50));
tags500 = tags500(tags500 <= 440);
prob500 = zeros(size(tags500));
for idx = 1:length(tags500)
    exp = sprintf("tags%dnumchannels500", tags500(idx));
    prob500(idx) = eval(exp);
end

% figure();
% hold on
% plot(tags50, prob50);
% plot(tags350, prob350);
% plot(tags500, prob500);
% legend(["50 Channels", "250 Channels", "500 Channels"]);
% title(["P(collision) for Varying Number" "of Frequency Channels"], FontSize=20);
% xlabel("Number of Tags", FontSize=18);
% ylabel("P(collision)", FontSize=18);
% set(gca, "YScale", "log")

figure();
hold on
plot(tags50 / 50, prob50, LineWidth=3);
plot(tags350 / 350, prob350, LineWidth=3);
plot(tags500 / 500, prob500, LineWidth=3);
legend(["50 Channels", "250 Channels", "500 Channels"]);
title(["Collisions for Varying Number" "of Frequency Channels"], FontSize=20);
xlabel("Activity Factor (\alpha)", FontSize=18);
ylabel("Collisions (%)", FontSize=18);

