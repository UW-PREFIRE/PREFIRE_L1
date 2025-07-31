%% PREFIRE Emission Precalculator
% running exitance subroutine for BRF is really slow, especially on SRF
% grid, since the Planck function is a fixed thing, this routine will
% calculate and store the 91 relevant temperatures for the calibration
% thermal_model_data
% 5/26/2020 B. added factor of pi to agree with A. blackbody calculations
% 5/31/2020 B. imported TIRS settings from modular file
% 11/13/2022 increase size of table to accommodate hot-ops up to 350K and
% increased resolution

[script_path, ~, ~] = fileparts(which('Emission_PreCalculator_v2'));
orig_path = path;
path(orig_path, fullfile(script_path, 'functions'));

pd = configure_toplevel_IO;

[TIRS] = prefire_TIRS(pd,2,0); %first argument is slit/pixel ratio, second is low/high gain (non-zero/zero)
channel_sampling = TIRS.d/TIRS.dispersion; %nominal pixel width divided by dispersion of a spectral micron per physical micron
delta_wave = channel_sampling / 100; % 6/700 microns
start_wave = delta_wave * 50;  % 0.4286 um
stop_wave = delta_wave * 7000; % 60.0 um
wavelen = (start_wave:delta_wave:stop_wave)'; 

T_space = 2.725;      % cold space scene temperature

% Reverting back to precomputing on a grid spacing of 1K. These
% values are linearly interpolating when applying to instrument
% data in L1B, so the 1K grid should be precise enough.
% Also, make the T-range wider to eventually accomodate lab test data.
start_temp = 100;
stop_temp = 350;
T_cal = start_temp:stop_temp;
num_temps = length(T_cal);

num_wave = size(wavelen, 1);
calibrator_emission = zeros([num_wave,num_temps]);
space_emission = zeros([num_wave,1]);

for i=1:num_wave
    %w2-w1 = 1.0 um gives emiss units W/m2/sr/um
    %w2-w1 = 0.1 um requires extra factor of 10 below
    w1 = abs(wavelen(i)-0.05); 
    w2 = wavelen(i)+0.05;
    % determine scene emission in wavelength range
    space_emission(i,1) =  10*10000*exitance(w1,w2,T_space)/pi;
    % factor of pi to agree with A. blackbody curves
    for n=1:num_temps 
        calibrator_emission(i,n) =  10*10000*exitance(w1,w2,T_cal(n))/pi;
    end
end %for i

%dumping cal_emission and space_emission means that total signal calculator will
%never need to call exitance.m
%% get the git hash from command line; this will be stored in the
%% MATFILE; we set the filename here with a manual version number &
%% date.

output_name = [ ...
    pd.instrument_model_dir, ...
    '/PREFIRE_Emission_preCalculated_v0.3.0_', ...
    datestr(datetime('now'), 'yyyy-mm-dd'), ...
    '.mat'];

[stat, git_short_hash] = system('git rev-parse --short HEAD');
% removes the trailing return character
git_short_hash = strip(git_short_hash);

save(output_name, ...
     'wavelen', 'T_space', 'T_cal', ...
     'calibrator_emission', 'space_emission', ...
     'git_short_hash');

path(orig_path);
