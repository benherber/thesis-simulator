%% CALCULATE

% 1−1/𝑁^𝐹*𝑁!/(𝑁−𝐹)!
clc; clear; close all;

slots = 10:20:110;
labels = strings(1, length(slots));

results = NaN(length(slots), slots(end));
for idx = 1:size(results, 1)
    num_slots = slots(idx);
    prob = @(family_size) 1 - (1 / (num_slots ^ family_size)) * ...
                (factorial(num_slots) / factorial(num_slots - family_size));
    families = 1:num_slots;
    results(idx, 1:length(families)) = arrayfun(prob, families);
    labels(idx) = sprintf("%d Slots", slots(idx));
end

%% PLOT

plot(1:slots(end), results);
xlim([0, 40]);
xlabel("Number of Tags in Family", FontSize=18);
ylabel("Probability of Collision", FontSize=18);
legend(labels);

%% CALCULATE

% 1−1/𝑁^𝐹*𝑁!/(𝑁−𝐹)!
clc; clear; close all;
f1 = figure;
f2 = figure;
ax1 = axes(f1);
ax2 = axes(f2);
hold(ax1, "on");
hold(ax2, "on");

slots = [8, 16, 32, 64];
labels = strings(1, length(slots));

for idx = 1:size(slots, 2)
    num_slots = slots(idx);
    prob = @(family_size) 1 - (1 / (num_slots ^ family_size)) * ...
                (factorial(num_slots) / factorial(num_slots - family_size));
    families = 1:(num_slots);
    result = arrayfun(prob, families);
    activities = families / num_slots;
    plot(ax1, activities, result);
    plot(ax2, families, result)
    labels(idx) = sprintf("%d Slots", slots(idx));
end

xlabel(ax1, "Activity Factor (\alpha)", FontSize=18);
ylabel(ax1, "Probability of Collision", FontSize=18);
legend(ax1, labels);
xlabel(ax2, "Number of Active Tags", FontSize=18);
ylabel(ax2, "Probability of Collision", FontSize=18);
legend(ax2, labels);
xlim([0, 30])
set(ax2, "YScale", "log")