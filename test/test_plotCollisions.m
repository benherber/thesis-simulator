clear; clc;
addpath(genpath(".."));
folder_name = "FinalSlurm";
figure();
xlabel("Number of Tags", FontSize=18);
ylabel("Probability of Collision", FontSize=18);
hold on

%% Aggregate PPersist
clearvars -except folder_name;

files = dir(fullfile("data/" + folder_name + "/", "ppersistcollisions*.mat"));
for idx = 1:length(files)
    load(fullfile(files(idx).folder, files(idx).name));
end

num_tags = ceil(linspace(1, 350, 50));
num_tags = num_tags(1:end-3);
probs = zeros(size(num_tags));
for idx = 1:length(probs)
    probs(idx) = eval(sprintf("tags%dppersist10", num_tags(idx)));
end

plot(num_tags, probs, LineWidth=3);

%% Aggregate Random Frequency
clearvars -except folder_name;

files = dir(fullfile("data/" + folder_name + "/", "randfreqhopcollisions*.mat"));
for idx = 1:length(files)
    load(fullfile(files(idx).folder, files(idx).name));
end

num_tags = ceil(linspace(1, 350, 50));
num_tags = num_tags(1:end-3);
probs = zeros(size(num_tags));
for idx = 1:length(probs)
    probs(idx) = eval(sprintf("tags%dnumchannels350", num_tags(idx)));
end

plot(num_tags, probs, LineWidth=3);

%% Aggregate Full Rand
clearvars -except folder_name;

files = dir(fullfile("data/" + folder_name + "/", "fullfreqppersistcollisions*.mat"));
for idx = 1:length(files)
    load(fullfile(files(idx).folder, files(idx).name));
end

num_tags = ceil(linspace(1, 350, 50));
num_tags = num_tags(1:end-3);
probs = zeros(size(num_tags));
for idx = 1:length(probs)
    probs(idx) = eval(sprintf("tags%dnumchannels350ppersist10", num_tags(idx)));
end

plot(num_tags, probs, LineWidth=3);

%% Plot specifics

legend(["P-Persistence", "Random FH", "Both"]);