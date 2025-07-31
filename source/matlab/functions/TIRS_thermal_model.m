function [thermistor_temps] = TIRS_thermal_model( ...
    pd, sc_latitude, sc_ctime, varargin)
%
% Thermal model for TIRS. The stored data are cubic spline fits
% to the temperature time series from the McKinely 2022 IEEE Aero
% conference paper (the "Cold Case" with more extreme temperature variation)
%
%
% input arguments:
% pd: structure from configure_toplevel_IO
% sc_latitude: column vector array with spacecraft latitude.
%     generally this would be a copy of Geometry/subsat_latitude
%     from a simulated L1 file.
% sc_ctime: continuous time array, same shape as sc_latitude.
%     generally this would be a copy of Geometry/ctime.
%
% optional input arguments:
% therm_noise_stdv: scalar specifying the random noise std. dev.
%     added to the modeled temperatures. In units of deg C.
%     Default is 0.02.
% rand_seed: seed number for random number generation, if noise is
%     added to the temperatures. The default uses a hardcoded
%     number, for reproducibility.

model_fpath = fullfile( ...
    pd.ancillary_data_dir, 'instrument_model', ...
    'TIRS_ColdCase_thermal_model_pp.h5');

if length(varargin) >= 1
    therm_noise_stdv = varargin{1};
else
    therm_noise_stdv = 0.02;
end
if length(varargin) == 2
    rand_seed = varargin{2};
else
    rand_seed = 177;
end


% Mckinley et al Thermal locations (Only a subset of those from the
% paper are included in the spline model - primarily because some
% of the temperature time series were almost identical. We have 7
% models for the 10 locations, 3 of them are very similar to the
% others.)

% 1 SC Interface
% 2 Lower Housing
% 3 Cal target
% 4 Toroid (Primary Mirror, Secondary mirror are similar)
% 5 Grating
% 6 Motor
% 7 Filter (Detector is similar)


% TIRS thermistor locations (these are the outputs). Order matters here.
% also, the mapping to Mckinley.
% 1 Toroidal Mirror -> 4 Toroid
% 2 Telescope       -> 2 Lower Housing (just a guess, not too
%                        important, not used in Inv. Inst. model)
% 3 filter          -> 7 Filter
% 4 cal target      -> 3 cal target
% 5 cal target      -> 3 cal target
% 6 cal target      -> 3 cal target
% 7 FPA center/back -> 7 Filter
% 8 FPA bracket     -> 7 Filter
%
% The McKinley data for Filter vs Detector was not discernable in
% the plot (the lines overlapped, with possibly 0.1 deg C
% difference.) We will add a small offset here (0.1) to make those
% two time series unequal.

% various (hardcoded, for now), temperature offsets
% again, we are assuming no large time gaps here (this acts on
% simulated data).
detector_temp_offset = 0.1;
% ad hoc temp adjustment to reduce the flux imbalance.
toroid_offset = 1.0;

num_frames = length(sc_latitude);

% some ad hoc stuff to add noise to the temps.
% This is highly temporally smoothed, and we scale the final result
% by the therm_noise_stdv value. The smoothing length is very ad
% hoc, has length scales of very roughly ~0.1-0.2 of the orbit;
% 1-sigma is ~500 samples of a ~8000 sample orbit period.
s = RandStream('mt19937ar', 'Seed', rand_seed);
therm_noise = randn(s, [8, num_frames]);
smoothing_kernel = exp( -(-5:0.002:5).^2 );
therm_noise_smooth = convn(therm_noise, smoothing_kernel, 'same');
mean_stdv = mean(std(therm_noise_smooth,0,2));
therm_noise_smooth = therm_noise_smooth * therm_noise_stdv/mean_stdv;

rel_ctime = sc_ctime - sc_ctime(1);
sc_phase = calc_orbit_normphase(sc_latitude, rel_ctime);
wrapped_sc_phase = mod(sc_phase, 1.0);

% toroid model: used in thermistor 1.
pp = mkpp( ...
    h5read(model_fpath, '/toroid/breaks'), ...
    h5read(model_fpath, '/toroid/coefs') );
toroid_temp = ppval(pp, wrapped_sc_phase);

% lower housing model: used in thermistor 2.
pp = mkpp( ...
    h5read(model_fpath, '/lower_housing/breaks'), ...
    h5read(model_fpath, '/lower_housing/coefs') );
lower_housing_temp = ppval(pp, wrapped_sc_phase);

% cal target model: used in thermistors 4,5,6
pp = mkpp( ...
    h5read(model_fpath, '/cal_target/breaks'), ...
    h5read(model_fpath, '/cal_target/coefs') );
cal_target_temp = ppval(pp, wrapped_sc_phase);

% filter model: used in thermistors 3,7,8
pp = mkpp( ...
    h5read(model_fpath, '/filter/breaks'), ...
    h5read(model_fpath, '/filter/coefs') );
filter_temp = ppval(pp, wrapped_sc_phase);

thermistor_temps = zeros([8, num_frames]);
thermistor_temps(1,:) = toroid_temp + toroid_offset;
thermistor_temps(2,:) = lower_housing_temp;
thermistor_temps(3,:) = filter_temp;
thermistor_temps(4,:) = cal_target_temp;
thermistor_temps(5,:) = cal_target_temp;
thermistor_temps(6,:) = cal_target_temp;
% this is used in InvInstrModel as detector temp; apply the
% small offset in order to produce non-zero irradiance
% for the related background terms.
thermistor_temps(7,:) = filter_temp + detector_temp_offset;
thermistor_temps(8,:) = filter_temp + detector_temp_offset/2;

% add noise
thermistor_temps = thermistor_temps + therm_noise_smooth;

end
