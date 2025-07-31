function [qflag] = create_L1_quality_flags(...
    calib_bitflags, detector_bitflags, channel_0_detector_bitflags, ...
    obs_bitflags);

% [qflag] = create_L1_quality_flags(...
%      calib_bitflags, detector_bitflags, channel_0_detector_bitflags,
%      obs_bitflags);
%
% Convert bit flags into int8 quality flags with 0,1,2 values.
% This function contains our qualitative assessment of the severity
% of each condition encoded by the bit flag states: in other words,
% which conditions create low quality or bad data.
%
% given input dimensions atrack, xtrack, spectral:
%
% Inputs (all assumed to be uint types)
%     calib_bitflags: (spectral, xtrack, atrack)
%     detector_bitflags: (spectral, xtrack)
%     channel_0_detector_bitflags: (xtrack,)
%     obs_bitflags: (atrack)
%
% outputs, all fields of a structure qflag, (all int8)
%     qflag.rad: (spectral, xtrack, atrack)
%     qflag.calib: (spectral, xtrack, atrack)
%     qflag.detector: (spectral, xtrack)
%     qflag.obs: (atrack)
%     qflag.channel_0_rad: (xtrack, atrack)
%     qflag.channel_0_detector (xtrack,)

[num_spec, num_xtrack, num_atrack] = size(calib_bitflags);

qflag = struct();
qflag.calib = zeros(num_spec, num_xtrack, num_atrack, 'int8');
qflag.detector = zeros(num_spec, num_xtrack, 'int8');
qflag.obs = zeros(num_atrack, 1, 'int8');
qflag.channel_0_rad = zeros(num_atrack, num_atrack, 1, 'int8');
qflag.channel_0_detector = zeros(num_xtrack, 1, 'int8');

% any conditions flagged in calib bitflags are qf2 (bad)
qflag.calib(calib_bitflags > 0) = 2;

% detector bit flags:
% mark 2, 3, 4, 5 as low quality, and
% only 0, 1 (masked or dead) as bad.
% Note bit 0 will be index 1 in MATLAB.
for bit = [3, 4, 5, 6]
    msk = bitget(detector_bitflags, bit) > 0;
    qflag.detector(msk) = 1;
    msk = bitget(channel_0_detector_bitflags, bit) > 0;
    qflag.channel_0_detector(msk) = 1;
end
for bit = [1, 2]
    msk = bitget(detector_bitflags, bit) > 0;
    qflag.detector(msk) = 2;
    msk = bitget(channel_0_detector_bitflags, bit) > 0;
    qflag.channel_0_detector(msk) = 2;
end

% Categorizing observation bitflags.  Most are in the low-quality category.
% Some are marked as bad:
%  [b2](3) large thermal/radiometric perturbation present (e.g., near eclipse
%       entrance)
%  [b5](6) large time interval between observation and the nearest calibration
%       sequence (due to corrupted or missing calibration sequence(s))
%  [b6](7) bus telemetry indicates that spacecraft attitude determination is
%       invalid
%  [b9](10) during a modeled-sun-avoidance type of bus slew
%  [b10](11) electronics warm-up period after powering instrument on
for bit = [1, 2, 4, 5, 8, 9]
    msk = bitget(obs_bitflags, bit) > 0;
    qflag.obs(msk) = 1;
end
for bit = [3, 6, 7, 10, 11]
    msk = bitget(obs_bitflags, bit) > 0;
    qflag.obs(msk) = 2;
end

% combine all
qflag.rad = zeros(num_spec, num_xtrack, num_atrack, 'int8');
qflag.channel_0_rad = zeros(num_xtrack, num_atrack, 'int8');

% matlab appears to only take 1D msk for indexing, so we have
% to loop through the xtrack columns of the detector flags.
for n = 1:num_xtrack
    qflag.rad(qflag.detector(:,n)==1,n,:) = 1;
    if qflag.channel_0_detector(n) == 1
        qflag.channel_0_rad(n,:) = 1;
    end
end
qflag.rad(:,:,qflag.obs==1) = 1;
qflag.channel_0_rad(:,qflag.obs==1) = 1;
% Not currently used, but left here for future changes in the calib
% bitflags. 
qflag.rad(qflag.calib==1) = 1;

for n = 1:num_xtrack
    qflag.rad(qflag.detector(:,n)==2,n,:) = 2;
    if qflag.channel_0_detector(n) == 2
        qflag.channel_0_rad(n,:) = 2;
    end
end
qflag.rad(:,:,qflag.obs==2) = 2;
qflag.channel_0_rad(:,qflag.obs==2) = 2;

qflag.rad(qflag.calib==2) = 2;

end



