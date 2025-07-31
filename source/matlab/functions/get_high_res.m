function [delta_wave, start_wave, stop_wave] = get_high_res(channel_sampling);
% this sets the high-res, subsampled wavelength grid to use for SRF
% calculations. Here we choose values that roughly divide each
% channel width (in v10 they are variable) by 100. The start and stop 
% should be beyond the edges of the first and last channels, in order 
% to get the tails correctly quantified
% switch to idealized data function if you want to see the impact
% of just the instrument function.
delta_wave = channel_sampling / 100; % 6/700 microns
start_wave = delta_wave * 50;  % 0.4286 um
stop_wave = delta_wave * 7000; % 60.0 um

end
