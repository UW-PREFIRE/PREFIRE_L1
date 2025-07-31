function [downlink_gaps] = TIRS_IIM_calc_downlink_gaps( ...
        ctime, downlink_ctime_origin, downlink_cadence, downlink_length)

% [downlink_gaps] = TIRS_IIM_calc_downlink_gaps(
%    ctime, downlink_ctime_origin, downlink_cadence,
%    downlink_length)
%
% compute downlink gaps that will occur in an input ctime range,
% given downlink parameters.
%
% Inputs:
% downlink_ctime_origin, downlink_cadence, downlink_length -
%    same inputs as TIRS_inverse_instrument_model
%
% returns:
% downlink_gaps : array of frame indices, with shape (n,2), where
%    n is the number of downlinks occuring in ctime.
%    Note that n could be zero, depending on the cadence and ctime
%    origin.

rel_ctime = ctime - downlink_ctime_origin;
first_downlink = ceil(rel_ctime(1) / downlink_cadence);
last_downlink = floor(rel_ctime(end) / downlink_cadence);
downlink_id = first_downlink:last_downlink;

num_downlinks = last_downlink - first_downlink + 1;
downlink_gaps = zeros([num_downlinks,2]);

for n=1:num_downlinks
    d = downlink_id(n);
    downlink_start = downlink_ctime_origin + d * downlink_cadence;
    downlink_stop =  downlink_start + downlink_length;
    idx = find((ctime >= downlink_start) & (ctime <= downlink_stop));
    downlink_gaps(n,1) = idx(1);
    downlink_gaps(n,2) = idx(end);
end

end
