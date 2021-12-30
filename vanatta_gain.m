% ----------------------------- %
% Find Gain of a Van Atta Array %
% Name   : vanatta_gain.m       %
% Author : Benjamin Herber      %
% Date   : Fall 2021            %
% ----------------------------- %

function resE = vanatta_gain(incoming_E, num_elements, fc, fb)
    tmp = 0;
    lambda_c = physconst("Lightspeed") / fc;
    lambda_b = physconst("Lightspeed") / fb;
    spacing = lambda_c / 2;

    for idx = 1:num_elements
        tmp = tmp + exp(1i * (find_phase(0, lambda_b, spacing, (idx - 1)) - find_phase(0, lambda_c, spacing, (idx - 1))));
    end

    resE = incoming_E * tmp;
end
