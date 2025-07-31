function [calseqs] = TIRS_IIM_prep_calseqs( ...
    ctime, calseq_ctime_origin, calseq_cadence, calseq_lengths )
% [calseqs] = TIRS_IIM_prep_calseqs(ctime, calseq_ctime_origin, ...
%             calseq_cadence, calseq_length)
%
% prepare MATLAB struct variable with data relating to calibration
% sequences (calseqs) in the TIRS data stream.
%
% Note that calseq are assumed to consist of:
%
% 1) transition from nadir (earth) view to space view
% 2) stare at space view
% 3) transition from space view port to internal cal view
% 4) stare at internal cal
% 5) transition from internal cal view back to nadir view
%
% Thus, there are 6 transition times, and 5 'dwell' periods in a
% single calseq.
%
% Inputs:
% ctime: a 'continuous time' array, generally /Geometry/ctime from
%     simulated L0 or L1 data file, of the frame integration times.
%     The ctime is assumed to have a few characteristics:
%     (Not necessarily by this function, but the following IIM
%     functions)
%     ctime separations (delta between successive samples) should
%     be larger than the assumed transition time (currently 0.1 s);
%     and smaller than the calseq_length input in this function.
%     (might need to be less than half the calseq_length - I did
%     not test this carefully, rather assuming that calseq length
%     is of order 4 or more ctime deltas.)
% calseq_ctime_origin: a scalar ctime value that serves as the
%     origin for the calseq timing. The actual value is not that
%     important; However, using the same ctime origin between a
%     series of simulated granules will ensure that the calseq
%     occurrences are evenly spaced across granule boundaries.
% calseq_cadence: number of seconds between calseq
% calseq_lengths: number of seconds in one of the calseq "stares" at
%     the space view and internal cal target. The total length of
%     the calseq is then:
%     calseq_lengths(1) + calseq_lengths(2) + 3 * calseq_transition_length
%     The transition length is current hard coded inside this
%     function at 0.1 s.
%
% outputs
% calseqs, a matlab struct array with fields:
%   event_lengths: time lengths of calseq events (5 elements)
%   event_relative_ctime: times of events relative to start of
%     calseq (6 elements)
%   time_length: total length in seconds of calseq
%   num: integer, number of calseqs within the span of ctime.
%   start_ctimes: array of num start ctimes, for each calseq.
%   end_ctimes: array of num end ctimes, for each calseq.
%   start_frames: array of num integer starting frame numbers
%   end_frames: array of num integer ending frame numbers.
%      start and end frame numbers are inclusive (e.g., matlab
%      convention, so that (start:end) includes all samples in the
%      calseq.
%   event_ctimes: array of (6, num) ctime values; this is
%      essentially the value of:
%      start_ctime + event_relative_ctime for each calseqs.


% hardcoded length of time of view of instrument internal housing
% (the transition time)
internal_look_length = 0.1;

calseqs.event_lengths = [ ...
    internal_look_length, calseq_lengths(1), ...
    internal_look_length, calseq_lengths(2), ...
    internal_look_length ];
calseqs.event_relative_ctime = [0, cumsum(calseqs.event_lengths)];
calseqs.time_length = calseqs.event_relative_ctime(end);

rel_ctime = ctime - calseq_ctime_origin;

% don't attempt to simulate a partial calseq at the end of ctime;
% so, subtract off the time length of one calseq from ctime(end).
first_calseq = ceil(rel_ctime(1) / calseq_cadence);
last_calseq = floor((rel_ctime(end) - calseqs.time_length) / calseq_cadence);

calseq_ct = first_calseq:last_calseq;
calseqs.num = length(calseq_ct);
calseqs.start_ctimes = calseq_ct * calseq_cadence + calseq_ctime_origin;
calseqs.end_ctimes = calseqs.start_ctimes + calseqs.time_length;

calseqs.start_frames = sorted_into_sorted(ctime, calseqs.start_ctimes);
calseqs.end_frames = sorted_into_sorted(ctime, calseqs.end_ctimes)+1;

calseqs.event_ctimes = zeros([6, calseqs.num]);
for n = 1:calseqs.num
    calseqs.event_ctimes(:,n) = calseqs.event_relative_ctime + ...
        calseqs.start_ctimes(n);
end

end


function [idx] = sorted_into_sorted(x, s)
% helper to mimic the operation of np.searchsorted.
% this version works for 1D vectors x and s, both assumed to be row
% vectors: shape is (1,n)). Assumes both x and s are sorted.
[~,sort_idx] = sort([x, s]);
num = length(x);
idx = find(sort_idx > num) - (1:length(s));
end
