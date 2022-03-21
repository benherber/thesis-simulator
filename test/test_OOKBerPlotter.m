clear; clc;
import bherber.thesis.Simulator bherber.thesis.TagType
addpath(genpath(".."));
filename = sprintf("data/ookbers::%s.mat", datestr(now, "mm-dd-yy-HH:MM:SS"));
save(filename);

errors = [];
try
    pool = gcp;
    num_threads = pool.NumWorkers;
    [params, num_symbs] = BerPlotterConstants;
    
    symbs_per_worker = (num_symbs + mod(num_symbs, num_threads)) / num_threads;
    
    tags = [1e-12;0;0];
    modes = [TagType.OOK];

    bers = NaN(1, num_threads);

    for snr_db = 0:10
        for thread = 1:num_threads
            scaled_params = BerPlotterConstants(symbs_per_worker);
            sim = Simulator(tags, modes, @(a) a, scaled_params, snr_db=snr_db);
            sim.auto_align()
            delay = sim.tag_delays(1);
            symbols = zeros(int32(scaled_params.simstep_sz), int32(scaled_params.sim_sym_ratio), 2);
            num_errors = 0;
    
            for jdx = 1:scaled_params.sim_sym_ratio
                    symbols(:, jdx, 2) = sim.step();
            end
    
            for idx = 1:symbs_per_worker
                symbols(:, :, 1) = symbols(:, :, 2);
                for jdx = 1:scaled_params.sim_sym_ratio
                    step = sim.step();
                    symbols(:, jdx, 2) = step;
                end
    
                curr_steps = symbols(:, :, 1);
                next_steps = symbols(:, :, 2);
                frame = [curr_steps(:); next_steps(:)];
                shifted_frame = delayseq(frame, -1 * delay);
                true_symbol = shifted_frame(1:scaled_params.symb_sz);
                demod_bit = sim.tags(1).demodulate(idx, true_symbol);
                if demod_bit ~= sim.tags(1).bits(idx)
                    num_errors = num_errors + 1;
                end
            end
    
            bers(thread) = num_errors / num_symbs;
            saveber = struct(sprintf("ber%dsnr%d", thread, snr_db), num_errors / num_symbs);
            save_out(filename, saveber);
        end
    
        savetotal = struct(sprintf("totalbersnr%d", snr_db), sum(bers, "omitnan"));
        save_out(filename, savetotal);
    end

catch errors
end

delete(gcp);
if ~isempty(errors)
    rethrow(errors);
end

%%
function save_out(filename, structure)
    save(filename, "-append", "-struct", "structure");
end
