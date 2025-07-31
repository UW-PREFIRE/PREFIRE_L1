function [wavelen, grating, filter, filter_ir, filter_er] = ...
                load_efficiency_data_v4(pd, start_wave, delta_wave, stop_wave);
% version 4 10/3/2022 bring in 'new' 1st order grating efficiency
%
% filter 4 data receieved 3/15/2021, in wavelength (microns), oversampled,
% no reflections, assume reflections are 30%
%
% new data received for filters 1,2,3 on February 15th, 2021, data is on a
% frequency grid, so it must be converted to microns for consistency with
% earlier code
%
% B. v2 returns 3 values in filter array for each interpolated wavelength
% the inner and outer reflectances are the first two values, the
% transmission is the third value
%
% function [wavelen, grating, filter] = load_efficiency_data( ...
%    start_wave, delta_wave, stop_wave);
%
% Load efficiency data (Filter transmissions from UReading,
%  grating efficiency data).
%
% Returns all data after using cubic interpolation to get the data
% all at a common wavelength grid. Note that the extrapolation
% value used in MATLAB interp1() is 0. This means that for a choice
% of start/stop_wave that goes outside the TIRS wavelength range
% will produce zero values in the filters and gratings for the
% out-of-range wavelengths.
%
% inputs:
% start_wave: start wavelength in um
% delta_wave: wavelength increment in um
% stop_wave: stop wavelength in um
%
% returns:
% wavelen: wavelength array, which will just be the MATLAB
%    result of start_wave:delta_wave:stop_wave
%    shape is (N,1)
% grating: grating efficiencies, shaped (N,5), for the different
%    orders contained in the file.
%    Index 3 (the -1 order) is the primary grating order used in
%    TIRS.
% filter: filter transmissions, shaped (N,5). Here some processing
%    is done to filter index 1, 5, following B. script.
%
%
% Note: initially, I think using a wavelength grid that evenly
% divides the TIRS channels is helpful. For the default parameters
% (as of Sept. 2019), these are:
% delta_wave = (54.0/63.0)/100;
% start_wave = delta_wave * 200; % 1.71 um
% stop_wave = delta_wave * 7000; % 60.0 um 


% Filter transmission mat file contains four arrays:
% filterf: wavelength grid (um) shape [4801,1]
% filter: filter response for 5 order sorting filters, shape
%   [5,4801]. Note that filter 5 appears to be a square
%   approximation?. Response given in pct.
% filter4f: wavelength grid (um) shape [696,1]
% filter4: filter response, shape [696,1]; this response is in
%   normalized units [0-1], but does go > 1 for the last few points
%   beyond 310 um wavelength
% B. code appears to replace filter(5,:) with the contents
% of filter4 and filter4f.

%load filter_3252019.mat

%on 1/23/2020, start replacing filter_3252019.mat with data received from
%G. Hawkins in PREFIRE_Spectra_4.xlsx, the data now includes both
%transmission and reflection as well as the diamond filter  these are
%spliced up into text files for reading

%on 2/15/2020, filter data for filters 1,2 and 3 were delivered on a
%different grids 
fpath_to_load = fullfile(pd.instrument_model_dir, 'filter_01222020f.txt');
filter0_wavelen = load(fpath_to_load, '-ascii'); %2-40 microns in nm units

fpath_to_load = fullfile(pd.instrument_model_dir, 'filter1_02152021f.txt');
filter1_wavenumber = load(fpath_to_load, '-ascii'); % wavenumber grid
filter1_wavelen = 10000000./filter1_wavenumber; %convert to nm

fpath_to_load = fullfile(pd.instrument_model_dir, 'filter23_02152021f.txt');
filter23_wavenumber = load(fpath_to_load, '-ascii'); % wavenumber grid
filter23_wavelen = 10000000./filter23_wavenumber; %convert to nm

fpath_to_load = fullfile(pd.instrument_model_dir, 'filter4_031521f.txt');
%filter4_wavelen = load('filter4_01222020f.txt', '-ascii'); %14-90 microns in nm units (truncated to match number of points in other files)
filter4_wavelen = load(fpath_to_load, '-ascii'); %oversampled data, in microns already

fpath_to_load = fullfile(pd.instrument_model_dir, 'filter0_01222020pct.txt');
filter0_dat = load(fpath_to_load, '-ascii'); %pct: inc. R, ex. R, T

fpath_to_load = fullfile(pd.instrument_model_dir, 'filter1_02152021pct.txt'); 
filter1_dat = load(fpath_to_load, '-ascii'); %pct: inc. R, ex. R, T 

fpath_to_load = fullfile(pd.instrument_model_dir, 'filter2_02152021pct.txt');
filter2_dat = load(fpath_to_load, '-ascii'); %pct: inc. R, ex. R, T

fpath_to_load = fullfile(pd.instrument_model_dir, 'filter3_02152021pct.txt');
filter3_dat = load(fpath_to_load, '-ascii'); %pct: inc. R, ex. R, T c

fpath_to_load = fullfile(pd.instrument_model_dir, 'filter4_031521.txt');
%filter4_dat = load('filter4_01222020pct.txt', '-ascii'); %pct: inc. R, ex. R, T
filter4_dat = load(fpath_to_load, '-ascii'); %oversampled data, only transmission, not percent

%load grating_12092018.mat %Grating efficiencies from Dan Wilson
%for order = -3,-2,-1(Main), 0 and 1
%
% grating_09082019f.txt contains the wavelength grid in nm, shaped
% (1596,1).
% grating_09082019g.txt contains the grating efficiencies for the
% different orders, shaped (1596,5).
fpath_to_load = fullfile(pd.instrument_model_dir, 'grating_09082019f.txt');
grating_wavelen = load(fpath_to_load, '-ascii');

fpath_to_load = fullfile(pd.instrument_model_dir, 'grating_09082019g.txt');
grating_dat = load(fpath_to_load, '-ascii');

% 1st order only grating efficiency
%fpath_to_load = fullfile(pd.instrument_model_dir, ...
%                               'PREFIRE Efficiency Prototype1 vs Design.xlsx');
%grating_1st_f = readmatrix(fpath_to_load, 'Sheet','Data', 'Range', 'A2:A267');
%grating_1st_g = readmatrix(fpath_to_load, 'Sheet','Data', 'Range', 'C2:C267');
%not assimilated yet...

% convert nm -> um
grating_wavelen = grating_wavelen / 1000.0;
filter0_wavelen = filter0_wavelen / 1000.0;
filter1_wavelen = filter1_wavelen / 1000.0;
filter23_wavelen = filter23_wavelen / 1000.0;
%filter4_wavelen = filter4_wavelen / 1000.0; %not needed with 3/15/21 data

% Interpolation to common grid:
wavelen = (start_wave:delta_wave:stop_wave)';
grating = zeros([size(wavelen,1), size(grating_dat,2)]);
grating_1st = zeros([size(wavelen,1),1]); %use same wavelength grid to interpolate new grating data
filter = zeros([size(wavelen,1), 5]);%size(filter_dat(:,1,1),2)]);
filter_ir = zeros([size(wavelen,1), 5]);%size(filter_dat(:,1,1),2)]);
filter_er = zeros([size(wavelen,1), 5]);%size(filter_dat(:,1,1),2)]);

for n = 1:4
    grating(:,n) = interp1(grating_wavelen, grating_dat(:,n), ...
                           wavelen, 'spline', 0.0);
end

%grating_1st(:,1) = interp1(grating_1st_f, grating_1st_g, wavelen, 'spline', 0.0); %use same wavelength grid to interpolate new grating data
grating(200:6350,3) = grating_1st(200:6350,1); %overwrite the 2-54 um range with newer 1st order grating info

% third dimension of the filter array is transmission
filter(:,1) = interp1(filter0_wavelen, filter0_dat(:,3), wavelen, 'spline', 0.0);
filter(:,2) = interp1(filter1_wavelen, filter1_dat(:,3), wavelen, 'spline', 0.0);
filter(:,3) = interp1(filter23_wavelen, filter2_dat(:,3), wavelen, 'spline', 0.0);
filter(:,4) = interp1(filter23_wavelen, filter3_dat(:,3), wavelen, 'spline', 0.0);
filter(:,5) = 100*interp1(filter4_wavelen, filter4_dat(:,1), wavelen, 'spline', 0.0);
% 1st dimension of the filter array is incidence reflection
filter_ir(:,1) = interp1(filter0_wavelen, filter0_dat(:,1), wavelen, 'spline', 0.0);                      
filter_ir(:,2) = interp1(filter1_wavelen, filter1_dat(:,1), wavelen, 'spline', 0.0);                      
filter_ir(:,3) = interp1(filter23_wavelen, filter2_dat(:,1), wavelen, 'spline', 0.0);                      
filter_ir(:,4) = interp1(filter23_wavelen, filter3_dat(:,1), wavelen, 'spline', 0.0);                      
filter_ir(:,5) = 30;
% 2nd dimension of the filter array is exitance reflection
filter_er(:,1) = interp1(filter0_wavelen, filter0_dat(:,2), wavelen, 'spline', 0.0);
filter_er(:,2) = interp1(filter1_wavelen, filter1_dat(:,2), wavelen, 'spline', 0.0);
filter_er(:,3) = interp1(filter23_wavelen, filter2_dat(:,2), wavelen, 'spline', 0.0);
filter_er(:,4) = interp1(filter23_wavelen, filter3_dat(:,2), wavelen, 'spline', 0.0);
filter_er(:,5) = 30;
% third dimension of the filter array is transmission

% convert pct -> normalized
filter = filter/100;
filter_ir = filter_ir/100;
filter_er = filter_er/100;
% AJM: should we do something for the negative filter values?

end
