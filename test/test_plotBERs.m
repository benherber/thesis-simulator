clear; clc;
addpath(genpath(".."));
% folder_name = "AdroitSimulations";
folder_name = "SlurmSimulations";

%% Aggregate OOK

files = dir(fullfile("data/" + folder_name + "/", "ookbers*.mat"));
for idx = 1:length(files)
    load(fullfile(files(idx).folder, files(idx).name));
end


snrs = linspace(0, 24, 13);
bers = zeros(3, length(snrs));
orders = [2, 4, 8];
for order = 1:length(orders)
    for ebno = 1:length(snrs)
        bers(order, ebno) = eval(sprintf("totalber%dMEbNo%d", orders(order), snrs(ebno)));
    end
end

plt(snrs, bers, orders, "M-PAM");

%% Aggregate FSK
clearvars -except folder_name;

files = dir(fullfile("data/" + folder_name + "/", "fskbers*.mat"));
for idx = 1:length(files)
    load(fullfile(files(idx).folder, files(idx).name));
end

snrs = linspace(0, 18, 10);
bers = zeros(3, length(snrs));
orders = [2, 4, 8];
for order = 1:length(orders)
    for ebno = 1:length(snrs)
        bers(order, ebno) = eval(sprintf("totalber%dMEbNo%d", orders(order), snrs(ebno)));
    end
end

plt(snrs, bers, orders, "M-FSK");

%% Plot

function plt(x, y, orders, type)
    figure();
    semilogy(x, y, LineWidth=2.5);

    labels = strings(1, length(orders));
    for idx = 1:length(labels)
        labels(idx) = sprintf("Order %d", orders(idx));
    end
    legend(labels);

    xlabel("Eb/No (dB)", FontSize=18);
    ylabel("Bit-Error Rate", FontSize=18);
    title(sprintf("Bit-Error of %s", type), FontSize=20);
    set(gca, "YGrid", "on");
end