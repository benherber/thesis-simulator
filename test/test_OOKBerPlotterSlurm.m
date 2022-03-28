clear; clc;
import bherber.thesis.Simulator bherber.thesis.TagType
addpath(genpath(".."));
filename = sprintf("data/ookbers::%s.mat", datestr(now, "mm-dd-yy-HH:MM:SS"));
save(filename);

errors = [];
try
    if isempty(gcp("nocreate"))
      pc = parcluster("local");
      pc.JobStorageLocation = "/tmp/";
      parpool(pc, 12)
    end
    pool = gcp;
    num_threads = pool.NumWorkers;
    modulation_order = uint16(str2double(getenv("SLURM_ARRAY_TASK_ID")));
    [params, num_symbs] = BerPlotterConstants(1, modulation_order);
    
    symbs_per_worker = ceil(num_symbs / num_threads);
    num_symbs = symbs_per_worker * num_threads;
    bits_per_symb = log2(params.m_ary_modulation);
    
    tags = [1;0;0];
    modes = [TagType.OOK];

    bers = NaN(1, num_threads);

    for snr_db = linspace(0, 24, 7)
        parfor thread = 1:num_threads
            scaled_params = BerPlotterConstants(symbs_per_worker, modulation_order);
            sim = Simulator(tags, modes, @(a) a, scaled_params, snr_db=snr_db, complex_noise=true);
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
                demod_bits = sim.tags(1).demodulate(idx, true_symbol);
                start = ((idx - 1) * bits_per_symb) + 1;
                stop = start + bits_per_symb - 1;
                expected = sim.tags(1).bits(start:stop).';
                num_errors = num_errors + sum(expected ~= demod_bits);
            end
    
            bers(thread) = num_errors / (num_symbs * bits_per_symb);
        end
%         frame_sz = 1000;
%         if mod(symbs_per_worker, frame_sz) ~= 0
%             error("invalid frame size.");
%         end
% 
%         parfor thread = 1:num_threads
%             scaled_params = BerPlotterConstants(symbs_per_worker, modulation_order);
%             sim = Simulator(tags, modes, @(a) a, scaled_params, snr_db=snr_db);
%             sim.auto_align()
%             delay = sim.tag_delays(1);
%             symbols = zeros(int32(scaled_params.simstep_sz), int32(scaled_params.sim_sym_ratio), frame_sz);
%             num_errors = 0;
%     
%             for jdx = 1:scaled_params.sim_sym_ratio
%                     symbols(:, jdx, frame_sz) = sim.step();
%             end
%     
%             for idx = 1:(symbs_per_worker / frame_sz)
%                 symbols(:, :, 1) = symbols(:, :, frame_sz);
%                 for kdx = 2:frame_sz
%                     for jdx = 1:scaled_params.sim_sym_ratio
%                         step = sim.step();
%                         symbols(:, jdx, kdx) = step;
%                     end
%                 end
%     
%                 for jdx = 1:(frame_sz - 1)
%                     curr_steps = symbols(:, :, jdx);
%                     next_steps = symbols(:, :, jdx + 1);
%                     frame = [curr_steps(:); next_steps(:)];
%                     shifted_frame = delayseq(frame, -1 * delay);
%                     true_symbol = shifted_frame(1:scaled_params.symb_sz);
%                     demod_bits = sim.tags(1).demodulate((((idx - 1) * frame_sz) + jdx), true_symbol);
%                     start = ((((idx - 1) * frame_sz) + jdx - 1) * bits_per_symb) + 1;
%                     stop = start + bits_per_symb - 1;
%                     expected = sim.tags(1).bits(start:stop).';
%                     num_errors = num_errors + sum(expected ~= demod_bits);
%                 end
%             end
%     
%             bers(thread) = num_errors / (num_symbs * bits_per_symb);
%         end
    
        savetotal = struct(replace(sprintf("totalber%dMEbNo%d", modulation_order, snr_db), ...
            "-", "NEG"), sum(bers, "omitnan"));
        save_out(filename, savetotal);
    end

catch errors
end

% delete(gcp);
if ~isempty(errors)
    rethrow(errors);
end

%%
function save_out(filename, structure)
    save(filename, "-append", "-struct", "structure");
end
