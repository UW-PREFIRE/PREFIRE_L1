function [encoder_pos] = TIRS_IIM_compute_calseq_encoder( ...
    ctime, calseqs, encoder_space, encoder_calib, encoder_nadir )
%
% [encoder_pos] = TIRS_IIM_compute_calseq_encoder( ...
%     ctime, calseqs, encoder_space, encoder_calib, encoder_nadir )
%
% generate time series of encoder positions during cal seq.
%
% inputs:
% ctime: array of continuous times for the frame readouts.
%     generally this is /Geometry/ctime from simulated L0/L1 data.
% calseqs: calseq structure (see TIRS_IIM_prep_calseqs).
% encoder_space, encoder_calib, encoder_nadir:
%     the scalar integer values, for the encoder positions at each
%     of those three ports.
%
% encoder_pos: array of same length as ctime, with integer (int16)
%     encoder position values.
%
% Basic algorithm: create time series of encoder positions at the
% edges of all transitions, where the mirror is slewing between
% fixed stares.
% then, linearly interpolate that (in time) back to the ctimes of
% the frame readouts. This results in 'intermediate' encoder values
% that might be sampled while the mirror is slewing between
% positions.


calseq_encoder_pos = repmat( ...
    [encoder_nadir, encoder_space, encoder_space, ...
     encoder_calib, encoder_calib, encoder_nadir], ...
    [1,calseqs.num]);
calseq_event_ctimes = reshape(calseqs.event_ctimes,[1,6*calseqs.num]);

% pad endpoints so we don't extrapolate. assume the system is in
% nadir pointing outside the first and last cal seq.
calseq_encoder_pos = [encoder_nadir, calseq_encoder_pos, encoder_nadir];
calseq_event_ctimes = [ctime(1), calseq_event_ctimes, ctime(end)];

encoder_pos = interp1(calseq_event_ctimes, calseq_encoder_pos, ctime);
encoder_pos = int16(round(encoder_pos));

end

                                                  