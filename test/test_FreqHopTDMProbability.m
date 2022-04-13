%% Model Activity Factor in FH-TDM Scenario

pattern_length = 6;
num_subbands = 6;
num_timeslots = 3;

activity_factor = linspace(0, 1, 46656);
probs = zeros(size(activity_factor));
for idx = 1:length(activity_factor)
    probs(idx) = collision_prob(pattern_length, num_subbands, num_timeslots, activity_factor(idx));
end
% probs = arrayfun(@(alpha) collision_prob(pattern_length, num_subbands, num_timeslots, alpha), activity_factor);

%% FUNCTIONS

function prob_collision = collision_prob(pattern_length, num_subbands, num_timeslots, activity_factor)
    arguments
        pattern_length (1, 1) double
        num_subbands (1, 1) double
        num_timeslots (1, 1) double
        activity_factor (1, 1) double {mustBeInRange(activity_factor, 0, 1, "inclusive")}
    end

    freq_family = num_subbands ^ pattern_length;
    prob_no_collisions = ones(1, (ceil(activity_factor * freq_family) - 1));
    for col = 1:(ceil(activity_factor * freq_family) - 1)
        prob_col_collide_numerator = (freq_family / num_subbands) - 1;
        prob_col_collide_denominator = num_timeslots * (freq_family - col);
        prob_col_collide = col * (prob_col_collide_numerator / prob_col_collide_denominator);
        prob_col_not_collide = min([0, 1 - prob_col_collide]);
        prob_no_collisions(col) = prob_col_not_collide;
    end
    
    prob_collision = 1 - (prod(prob_no_collisions) ^ pattern_length);
end