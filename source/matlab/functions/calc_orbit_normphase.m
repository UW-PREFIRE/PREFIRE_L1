function [sc_phase] = calc_orbit_normphase(sc_latitude, relative_ctime)
%
% [sc_phase] = calc_orbit_normphase(sc_latitude, relative_ctime)
%
% convert spacecraft latitude (the ground track latitude) to
% normalized orbit phase, defined in the range [0,1].
% norm orbit phase = 0 as the spacecraft leaves the ascending node,
% and 1 as the spacecraft arrives at the ascending node.
%
% the method used here is a bit ad hoc and approximate, as we only
% need this for the instrument thermal model; do not use this
% function for any purpose were accuracy is needed (e.g., anything
% related to orbit position or geolocation)
%
% Since this is assumed to work on simulated data, the method is
% not robust to time gaps.
%
% inputs:
% sc_latitude: column vector array with spacecraft latitude
% relative_ctime: continuous time array, same shape as sc_latitude.
%     this should be relative (e.g., relative to the first point),
%     just to keep the magnitude of the values lower.
%
% returns:
% sc_phase: column vector, same shape as input arrays, with
%     normalized orbit phase.

dlat = diff(sc_latitude);
ascending_msk = dlat > 0;

% we need to find at least two points: a node, or either extrema in
% latitude.
anode_msk = (sc_latitude(2:end) > 0) & (sc_latitude(1:end-1) < 0);
dnode_msk = (sc_latitude(2:end) < 0) & (sc_latitude(1:end-1) > 0);

latmin_msk = (dlat(2:end) > 0) & (dlat(1:end-1) < 0);
latmax_msk = (dlat(2:end) < 0) & (dlat(1:end-1) > 0);

% this is not strictly correct (orbit eccentricity, earth ellipsoid
% shape, would alter these slightly, but this should be close enough.)
anode_phase = 0.0;
latmax_phase = 0.25;
dnode_phase = 0.5;
latmin_phase = 0.75;

% collect times and phases of these four events (asc/desc node, and
% latitude min/max.)
[t_tiepts_1, p_tiepts_1] = find_timephase_tiepts( ...
    relative_ctime, anode_msk, anode_phase);

[t_tiepts_2, p_tiepts_2] = find_timephase_tiepts( ...
    relative_ctime, latmax_msk, latmax_phase);

[t_tiepts_3, p_tiepts_3] = find_timephase_tiepts( ...
    relative_ctime, dnode_msk, dnode_phase);

[t_tiepts_4, p_tiepts_4] = find_timephase_tiepts( ...
    relative_ctime, latmin_msk, latmin_phase);

t_tiepts = [t_tiepts_1; t_tiepts_2; t_tiepts_3; t_tiepts_4];
p_tiepts = [p_tiepts_1; p_tiepts_2; p_tiepts_3; p_tiepts_4];

if length(t_tiepts) < 2
    error('Need at least 2 tiepoints in spacecraft latitude');
end

[t_tiepts, sorter] = sort(t_tiepts);
p_tiepts = p_tiepts(sorter);

% unwrap phase
p_tiepts = unwrap_phase(p_tiepts);

% linear fit/interpolation back to input ctime
linear_coef = polyfit(t_tiepts, p_tiepts, 1);
sc_phase = linear_coef(1) * relative_ctime + linear_coef(2);

end


function [t_tiepts, p_tiepts] = find_timephase_tiepts(ctime, msk, p_value);
% helper to extract 'tie points' in time and phase, given an input
% mask array that marks the locations of the tie points. The mask
% values are assumed to mark the time range specified by time
% values (i) and (i+1), hence we average those two values to get
% the tiepoint time. phase values at the tie points are constant,
% equal to the input p_value.
if any(msk);
    idx = find(msk);
    npts = length(idx);
    t_tiepts = 0.5*(ctime(idx)+ctime(idx+1));
    p_tiepts = ones(npts,1)*p_value;
else
    t_tiepts = [];
    p_tiepts = [];
end
end

function [phase_uw] = unwrap_phase(phase)
% helper function to unwrap the normalized phase
% cant use the builtin unwrap, unfortunately, as it assumes
% radians, and also attempts to unwrap both positive and negative
% shifts (whereas here there are only negative drops, because phase
% is alwaus monotonically increasing.)
p_diff = diff(phase);
phase_uw = phase;
wrap_idx = find(p_diff < 0);
for n=1:length(wrap_idx)
    phase_uw(wrap_idx(n)+1:end) = phase_uw(wrap_idx(n)+1:end) + 1;
end
end
