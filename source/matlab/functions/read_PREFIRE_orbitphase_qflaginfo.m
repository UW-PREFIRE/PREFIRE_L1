function [error_ca, qflag_info] = ...
        read_PREFIRE_orbitphase_qflaginfo(data_ctime, info_file)
% [error_ca, qflag_info] = ...
%     read_PREFIRE_orbitphase_qflaginfo(data_ctime, info_file)
%
% read information from csv file describing the sections of the
% orbit that should be flagged low quality or bad due to eclipse
% phenomena, or another other features that happen at regular orbit
% positions (phases) that drift seasonally.
% Currently this is hardcoded to define 4 orbit phase intervals.
% This function handles the temporal interpolation between two
% orbit phase range definitions. For example, if:
%
% ctime, level, orbit phase range
% 1000,      2, 100, 200
% 2000,      2, 120, 220
%
% and the requested data_ctime is 1500, then the returned interval
% will be 110 - 210 in orbit phase, with level 2.
% the returned level is the minimum of the two "bracketing"
% entries. In this way, changing the level to a zero, would
% effectively disable that orbitphase range from modifying the
% quality flag. For example,
%
% ctime, level, orbit phase range
% 1000,      1, 100, 200
% 1500,      1, 100, 200
% 2000,      0, 100, 200
%
% any ctime request between 1000 and 1500 will produce an interval
% 100-200 with level 1; and any ctime request between 1500 and 2000
% will be 100-200 with level 0.
%
% Inputs:
%     data_ctime: scalar ctime value for the requested orbit.
%         assumed to be the start of the orbit (granule_beg_ts).
%     info_file: string filepath to info CSV file.
%
% Outputs:
%     error_ca: error cell array
%     qflag_info: struct with following fields:
%        levels: 1D array of severity levels (0, 1 or 2), 4 values.
%          values. Level 0 effectively means it will be ignored.
%        orbit_phase_ranges: 2D array with shape (4,2), with the
%          orbit phase intervals, linearly interpolated from the
%          stored CSV data table values.

error_ca = {'#NONE#', 0};  % Default

% Note this reads into a custom MATLAB Table object
T = readtable(info_file);

qflag_info.level = zeros([4,1]);
qflag_info.ref = zeros([4,1]);
qflag_info.orbitphase_offsets = zeros([4,2]);

% locate time position inside info file
% throw error if the data ctime is larger than all times in info
% file, meaning the requested data has moved beyond the end of the table.
idx = find(T.ctime > data_ctime, 1, 'first');
if isempty(idx)
    error_ca = {['No data in orbitphase qflag is defined for this ' ...
                 'granule time'], 1};
    return
end

if (idx == 1)
    wts = [1., 0.];
    idx_b = 1;
    idx_e = 2;
else
    ctime_duration = T.ctime(idx) - T.ctime(idx-1);
    wt = (T.ctime(idx) - data_ctime) / ctime_duration;
    wts = [wt, 1-wt];
    idx_b = idx-1;
    idx_e = idx;
end

for n = 1:4
    c = 3+(n-1)*5;

    set_ID = int2str(n);

    col_name = ['qlevel_' set_ID];
    qflag_info.level(n) = T.(col_name)(idx_b);

    col_name = ['ref_' set_ID];
    ref = T.(col_name)(idx_b);
    if strcmp(ref, 'ref_to_eclipse_entrance')
       qflag_info.ref(n) = 11;
    elseif strcmp(ref, 'ref_to_eclipse_exit')
       qflag_info.ref(n) = 12;
    else
       qflag_info.ref(n) = 10;
    end

    % Note that the application of wts is a dot product on purpose
    for ic=1:2
       c_offset = c+1+ic;
       qflag_info.orbitphase_offsets(n,ic) = ...
                                     wts * table2array(T(idx_b:idx_e,c_offset));
    end

end

end
