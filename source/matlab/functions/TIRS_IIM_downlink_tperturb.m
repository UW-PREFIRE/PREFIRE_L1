function [thermistor_temps] = TIRS_IIM_downlink_tperturb( ...
    ctime, thermistor_temps, downlink_gaps, ...
    downlink_tperturb_size, downlink_tperturb_length);

[num_downlinks, ~] = size(downlink_gaps);
num_frames = length(ctime);
% note that 5 e-folds gets to 0.0067, and the resolution of
% the thermistors is about 0.02. So, 5 e-folds should be smaller
% than the resolution assuming the max T perturbation (right at the
% end of the downlink gap) is only ~1-2 K.
num_efolds = 5;

for n = 1:num_downlinks
    if downlink_gaps(n,2) < num_frames
        % exponentially asymptoting tperturb during downlink
        rel_ctime = ctime - ctime(downlink_gaps(n,1));
        start_fr = downlink_gaps(n,1);
        end_fr = downlink_gaps(n,2);
        rel_ctime_chunk = rel_ctime(start_fr:end_fr) - rel_ctime(start_fr);
        tperturb1 = downlink_tperturb_size * ...
            (1 - exp(-rel_ctime_chunk/downlink_tperturb_length));
        % Note that we apply temp pertubation to toroid only (slot 1)
        thermistor_temps(1,start_fr:end_fr) = ...
            thermistor_temps(1,start_fr:end_fr) + tperturb1';

        % exponentially decaying tperturb at end of downlink.
        rel_ctime = ctime - ctime(downlink_gaps(n,2));
        start_fr = downlink_gaps(n,2) + 1;
        end_fr = start_fr;
        end_time = num_efolds * downlink_tperturb_length;
        while (rel_ctime(end_fr) <= end_time) & (end_fr < num_frames)
            end_fr = end_fr + 1;
        end
        rel_ctime_chunk = rel_ctime(start_fr:end_fr) - rel_ctime(start_fr);
        tperturb2 = exp(-rel_ctime_chunk/downlink_tperturb_length) * ...
            tperturb1(end);
        % Note that we apply temp pertubation to toroid only (slot 1)
        thermistor_temps(1,start_fr:end_fr) = ...
            thermistor_temps(1,start_fr:end_fr) + tperturb2';

    end
end

end
