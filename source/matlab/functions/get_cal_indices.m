function [cal_indices] = get_cal_indices(seconds,cadence);
%this function determines the scenes within a simulated data run that will
%be replaced with simulated calibration data to mimic TIRS operation
offset = 1;
num_cals = floor(size(seconds,2)/cadence);
cal_indices = (0:num_cals)*cadence+offset;

end
