classdef HopsetsGenerator
    %HOPSET_GENERATOR enabling efficient frequency hopping.
    %   Gernerate optimal hopsets for Multi-Tag Communication.

    methods (Static)
        function hopsets = generate(num_sets, num_channels, num_conflicts_allowed)
            %GENERATE a specified number of optimal hopsets.
            %   If more hopsets are requested than channels, conflicts must
            %   be allowed. Conflicts will determine maximum number of
            %   allowed sets. Shorter pattern length will be favored.
            %   
            %   Alorithm Invarients:
            %     1. Unless repeated element, columns will be ascending
            %     from top to bottom
            %     2. Repeats occur in a diagonal fashion

            arguments
                num_sets double {mustBePositive}
                num_channels double {mustBePositive}
                num_conflicts_allowed double ...
                    {mustBeInRange(num_conflicts_allowed, 0, 2, "inclusive")} = 0;
            end

            if (num_conflicts_allowed == 0) || (num_sets <= num_channels)
                if num_sets > num_channels; error("Cannot find hopsets w/o conflicts"); end
                hopsets = zeros(num_sets, num_channels);
                base = 1:num_channels;
                for idx = 1:num_sets
                    hopsets(idx, :) = base(1:num_channels);
                    base = [base(end), base(1:(num_channels - 1))];
                end
            elseif num_conflicts_allowed == 1
                if num_channels < 5; error("To few useable channels."); end
                optimal_k = 3; % FIXME: Favor smaller k
%                 base = zeros(2 * optimal_k, optimal_k);
                base = [ ...
                        [1, 3, 2]; ...
                        [1, 4, 5]; ...
                        [3, 2, 4]; ...
                        [4, 2, 1]; ...
                        [2, 5, 3]; ...
                        [5, 1, 3]; ...
                       ];
%                 for idx = 1:(2 * optimal_k)
%                     base(idx, :) = base(idx, :) + idx;
%                 end
%                 for idx = 1:optimal_k
%                     start = (2 * idx);
%                     stop = (2 * optimal_k);
%                     adjusted = base(:, idx) - 1;
%                     base(start:stop, idx) = adjusted(start:stop);
%                 end
%                 for idx = 1:optimal_k
%                     if optimal_k == 2
%                         base((2 * idx) - 1, idx) = 1;
%                         base(2 * idx, idx) = 1;
%                     else
%                         base((2 * idx) - 1, idx) = idx;
%                         base(2 * idx, idx) = idx;
%                     end
%                 end
%                 for idx = 1:(2 * optimal_k)
%                     for tdx = 1:optimal_k
%                         if base(idx, tdx) ~= 0; continue; end
%                         base(idx, tdx) = find_lowest(base(idx, :), base(:, tdx), num_channels); 
% 
%                         if isnan(base(idx, tdx))
%                             base(idx, tdx) = base(idx - 1, tdx);
%                             base(idx - 1, tdx) = 0;
%                             base(idx - 1, tdx) = find_lowest(base(idx - 1, :), base(:, tdx), num_channels);
%                         
%                             if isnan(base(idx - 1, tdx))
%                                 base(idx - 1, tdx) = find_lowest([], base(:, tdx), num_channels);
%                             end
%                         end
%                     end
%                 end

                need_len = floor(num_sets / (2 * optimal_k));
                hopsets = zeros(need_len * (2 * optimal_k), optimal_k);
                for idx = 1:need_len
                    channels_needed = (2 * optimal_k) - 1;
                    start = (idx - 1) * (2 * optimal_k) + 1;
                    stop = idx * (2 * optimal_k);
                    hopsets(start:stop, :) = base + ((idx - 1) * channels_needed);
                end
                channels_used = need_len * (2 * optimal_k - 1);
                channels_left = num_channels - channels_used;
                sets_left = num_sets - need_len * (2 * optimal_k);
                if sets_left ~= 0
                    if sets_left > channels_left; error("Too few available subbands."); end
                    final_hopsets = zeros(sets_left, optimal_k);
                    base = (channels_used + 1):(channels_used + channels_left);
                    for idx = 1:sets_left
                        final_hopsets(idx, :) = repelem(base(idx), optimal_k);
                    end
                    hopsets = [hopsets; final_hopsets];
                end
            end

            function res = find_lowest(a, b, hi)
                res = NaN;
                for jdx = 1:hi
                    in_a = sum(ismember(a, jdx));
                    in_b = sum(ismember(b, jdx));
                    if ~(in_a || in_b)
                        res = jdx;
                        break;
                    end
                end
            end
        end
    end
end