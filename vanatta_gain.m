% % ----------------------------- %
% % Find Gain of a Van Atta Array %
% % Name   : vanatta_gain.m       %
% % Author : Benjamin Herber      %
% % Date   : Fall 2021            %
% % ----------------------------- %
% 
% function res = vanatta_gain(num_elements, fc, fb)
% %VANATTA_GAIN of a backscatter tag in a Van Atta Configuration.
% %   Calculate gain of a backscatter tag in a Van Atta Configuration for
% %   given input and output frequencies.
% 
%     tmp = 0;
%     lambda_c = physconst("Lightspeed") / fc;
%     lambda_b = physconst("Lightspeed") / fb;
%     spacing = lambda_c / 2;
% 
%     for idx = 1:num_elements
%         tmp = tmp + exp(1i * (find_phase(0, lambda_b, spacing, (idx - 1)) ...
%             - find_phase(0, lambda_c, spacing, (idx - 1))));
%     end
% 
%     res = tmp;
% end
