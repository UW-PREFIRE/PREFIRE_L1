function [TIRS] = prefire_TIRS(ROIC_tau, pd, varargin)
if(length(varargin)>0)
    slit_flag=varargin{1};
else
    slit_flag=2; %default value
end

if(length(varargin)>1) %B. removed dependence of noise on input choice 10/2/2022, leaving in flag for reverse compatibility
    noise_flag=varargin{2};
else
    noise_flag = 1; %'high_gain' and low noise is default
end

if(length(varargin)>2)
    TIRS_flag=varargin{3};
else
    TIRS_flag = 0; % zero indicates model parameters, 1 and 2 refer to real TIRS measured values for instruments 1 and 2
end
% the nature of the interaction between varargin and pd means that
% we can easily invoke prefire_TIRS incorrectly, and get no errors,
% and wind up with incorrect data.
% Use the assertion to make sure pd is actually a struct.
assert(isstruct(pd));

% these are focal plane and ROIC design parameters
TIRS.spectral_channels = 64;
TIRS.spatial_scenes = 8;
TIRS.d = 180; %(um) width and length of detector 
TIRS.d2 = 370;%365; % pixel strip gap in microns changed from 365 on 7/12/2020 
TIRS.tau = ROIC_tau; % integration time in seconds
TIRS.cadence = 180; % number of images between starts of calibration sequences
% encoder position ( 0,26,46,66 )*16384 / 360 + offset (0 is safe, 20 is nadir, 40 is internal target, 60 is space)
TIRS.nadir_enc = 3148; % encoder position for nadir pointing
TIRS.space_enc = 4866; % encoder position for space pointing
TIRS.cal_enc = 3962;  % encoder position for calibration target pointing

% per B. Ho 6/3/2022, the final gain is 20.2e6 counts/volt because R_SCLK_DIV = 1 (or divide by 4) giving chopper frequency = 62.5 kHz (mclk = 125 kHz), and sliding window = 1, although there will be some gain variation.
% "As the gain is proportional to integration time, chopper frequency, and
% 2^(- sliding window), it becomes 20.5e6 counts/V for integration time = 694 ms, chopper frequency = 62.5 kHz and sliding window = 1, which I think is the intended flight configuration."

% B. will remove dependence on 'noise_flag' and fix the gain 10/2/2022
TIRS.g = 20.2e6; % alternating sign will be accounted for in ROIC_map
% B. not sure what to set for noise, guess here
TIRS.noise = 90e-9/sqrt(TIRS.tau);
%if (noise_flag ~= 0)
    %TIRS.noise = 60e-9/sqrt(TIRS.tau); % 60 nV/Hz^1/2 noise level for high gain
%    TIRS.noise = 120e-9/sqrt(TIRS.tau); % latest CBE 2/25/2021
%    TIRS.g = -7.6e6; %(counts/Volt) ROIC gain for high gain
%else 
    %TIRS.noise = 90e-9/sqrt(TIRS.tau); % 90 nV/Hz^1/2 noise level for low gain
%    TIRS.noise = 180e-9/sqrt(TIRS.tau); % guess at new CBE 2/25/2021
%    TIRS.g=-7.6e6/3; %(counts/Volt) ROIC gain for low gain
%end
TIRS.o = 24000; %offset counts

%mask for L1 algorithm to ignore disconnected/bad pixels
%initialize with non-blocked/blocked channels common to all scenes
mask = ...
[  1  1  0  0  1  1  1  1  0  0  1  1  1  1  1  1  1  0  0  1];
%  0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19
mask = [mask, ...
   1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  0  0  1  1  1];
% 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39
mask = [mask, ...
   1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1];
% 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59
mask = [mask, 1, 1, 1, 1];
TIRS.mask_blocked = [mask; mask; mask; mask; mask; mask; mask; mask];

if (TIRS_flag == 1)
% measurements for each pixel provided as D* values from blackbody testing
% for focal plane 001, SN101, data provided in 
% "FPA101 Test Results 2022-04-17.xlsx"
    TIRS.fl = 8100; % focal length of optics in microns, spatial dimension
    TIRS.fl2 = 6471; % focal length of optics in microns, spectral dimension
    fpath_to_load = fullfile(pd.instrument_model_dir, ...
                             'FPA101 Test Results 2022-04-17.xlsx');
    TIRS.Dstar = readmatrix(fpath_to_load, 'Sheet', 'Pixel Map', 'Range', ...
                            'D3:K66');
    TIRS.Dstar = 1e8*TIRS.Dstar; %table was in Jones / 10^8
    TIRS.R = TIRS.noise.*TIRS.Dstar./(1e-4*TIRS.d*sqrt(TIRS.tau));
    TIRS.x_fp = -90.0; % offset of focal plane in spectral dimension for instrument 1
    TIRS.tilt_fp = deg2rad(0.0); %tilt of focal plane, rotation about x-axis, nominal value determined acceptable after alternate explanation of test data (B. June 2023)
    TIRS.y_fp = 0.0;%50.0; % offset of focal plane in spatial dimension
%first argument is slit/pixel ratio, 
%second is low/high gain (zero/non-zero) [defunct]
%third is instrument # (zero or no entry is model)
    %mask out bad pixels
    TIRS.mask = TIRS.mask_blocked;
    TIRS.mask(2,13) = 0;TIRS.mask(2,39)=0;TIRS.mask(2,54)=0;TIRS.mask(2,64) = 0;
    TIRS.mask(3,13:14)=0;TIRS.mask(3,21:23)=0;TIRS.mask(3,44) = 0; TIRS.mask(3,62) = 0;
    TIRS.mask(5,31)=0;
    TIRS.mask(6,5)=0;TIRS.mask(6,29)=0;TIRS.mask(6,52:53)=0;TIRS.mask(6,55)=0;
    TIRS.mask(7,22)=0;TIRS.mask(7,24:26)=0;TIRS.mask(7,39)=0;TIRS.mask(7,58)=0;TIRS.mask(7,61) = 0; 
    TIRS.mask(8,7:8)=0;TIRS.mask(8,12)=0;TIRS.mask(8,23)=0;TIRS.mask(8,28)=0;TIRS.mask(8,58)=0;
    %empirically adjusted values to scale signal
    arb = [1,1,1,1,13.6820809547881,3.82135567376512,2.29636244482096,1.24212089023718,1,1,1.03322424833353,1.02273276275862, ones(64-12,1)'];
    TIRS.arb = [arb; arb; arb; arb; arb; arb; arb; arb]';
    TIRS.centers = [0,0,1.30000000000000,2.10000000000000,7.70000000000000,6,6,6,6.33000000000000,7.17000000000000,9,9.60000000000000,9.83000000000000,10.9000000000000,11.4000000000000,12.2000000000000,13,13.9000000000000,14.8000000000000,15.8000000000000,16.4500000000000,17.3000000000000,18.1400000000000,18.9900000000000,19.8300000000000,20.6700000000000,21.5200000000000,22.3600000000000,23.2000000000000,24.0500000000000,24.8900000000000,25.7400000000000,26.5800000000000,27.4200000000000,28.2700000000000,29.1100000000000,29.9600000000000,30.8000000000000,31.6400000000000,32.4900000000000,33.3300000000000,34.1700000000000,35.0200000000000,35.8600000000000,36.7100000000000,37.5500000000000,38.3900000000000,39.2400000000000,40.0800000000000,40.9200000000000,41.7700000000000,42.6100000000000,43.4600000000000,44.3000000000000,45.1400000000000,45.9900000000000,46.8300000000000,47.6800000000000,48.5200000000000,49.3600000000000,50.2100000000000,51.0500000000000,51.8900000000000,52.7400000000000];
    % central wavelength values from instrument model, currently shifted 0.4 um down from nominal

elseif (TIRS_flag == 2)
% measurements for each pixel provided as D* values from blackbody testing
% for focal plane 002, SN102, data provided in 
% "FPA102 Test Results.xlsx"
    TIRS.fl = 8100; % focal length of optics in microns, spatial dimension
    TIRS.fl2 = 6471; % focal length of optics in microns, spectral dimension
    fpath_to_load = fullfile(pd.instrument_model_dir, ...
                             'FPA102 Test Results.xlsx');
    TIRS.Dstar = readmatrix(fpath_to_load, 'Sheet', 'Pixel Map', 'Range', ...
                            'D3:K66');
    TIRS.Dstar = 1e8*TIRS.Dstar; %table was in Jones / 10^8
    TIRS.R = TIRS.noise.*TIRS.Dstar./(1e-4*TIRS.d*sqrt(TIRS.tau));
    
    %there is a +0.3 spectral micron shift evident in the spectral
    %calibration data, this is (0.3/0.84)*180 physical microns 
    %HOWEVER, the field points (in all spatial spectral cal) were inadvertently shifted
    %in elevation by about 1.0 degrees, this introduces an additional
    %spectral shift of (1.0/theta_spectral)*180 physical microns
    %in cross track 180 um = 1.3 degrees, scale by fno 180 um = 1.3*2/1.6
    %=1.63 degrees so 1.0 degree elevation shift is 180(1/1.6) = 112.5
    %microns!!  the total shift is 112.5+64 = 176.5 microns
    TIRS.x_fp = 112.5+64.0 ;%   
    TIRS.tilt_fp = deg2rad(0.0);%1.0; %tilt of focal plane, rotation about x-axis
    TIRS.y_fp = 0.0;%50.0; % offset of focal plane in spatial dimension
    %mask out bad pixels
    %Bad pixels from resistance measurement: A 39, 64; C 33-35, 38-64; E 40, 44, 47, 48, 53, 55, 58, 60, 62-64; G 38, 41, 48, 60; H 41, 51, 53, 58, 60, 61
    TIRS.mask = TIRS.mask_blocked;
    TIRS.mask(1,39) = 0;TIRS.mask(1,64) = 0;
    TIRS.mask(3,33:35) = 0; TIRS.mask(3,38:64) = 0;
    TIRS.mask(5,40) = 0;TIRS.mask(5,44) = 0;TIRS.mask(5,47) = 0;TIRS.mask(5,48) = 0;TIRS.mask(5,53) = 0;TIRS.mask(5,55) = 0;TIRS.mask(5,58) = 0;TIRS.mask(5,60) = 0;TIRS.mask(5,62:64) = 0;
    TIRS.mask(7,38) = 0;TIRS.mask(7,41) = 0;TIRS.mask(7,48) = 0;TIRS.mask(7,60) = 0;
    TIRS.mask(8,41) = 0;TIRS.mask(8,51) = 0;TIRS.mask(8,53) = 0;TIRS.mask(8,58) = 0;TIRS.mask(8,60) = 0;TIRS.mask(8,61) = 0;
    %additional troublesome pixels identified in radcal
    TIRS.mask(1,48) = 0;
    TIRS.mask(4,52) = 0;
    TIRS.mask(5,45)=0;TIRS.mask(5,50) = 0;
    TIRS.mask(6,39)=0;TIRS.mask(6,44)=0;TIRS.mask(6,52)=0;
    TIRS.mask(7,56)=0;TIRS.mask(7,57)=0;TIRS.mask(7,62)=0;
    arb = [1,1,1,1,6.17875663999124,1.73409602098750,1.14788776468947,1.11549047981300,1,1,0.981211191069921,0.979230989401896,ones(64-12,1)'];
    TIRS.arb = [arb; arb; arb; arb; arb; arb; arb; arb]';
    %these centers were taken from wl structure out of geometry_v3 and then
    %manually adjusted with radcal data in radcal.invert.m
    TIRS.centers = [0.827402116402116,1.67121164021164,2.51502116402116,3.35883068783069,7.70264021164021,5.99644973544974,6.24025925925926,6.73406878306878,7.57787830687831,8.52168783068783,9.36549735449735,10.1093068783069,10.9531164021164,11.7469259259259,12.5907354497354,13.4845449735450,14.2783544973545,15.1721640211640,16.0159735449735,16.8597830687831,17.7035925925926,18.5474021164021,19.3912116402116,20.2350211640212,21.0788306878307,21.9226402116402,22.7664497354497,23.6102592592593,24.4540687830688,25.2978783068783,26.1416878306878,26.9854973544974,27.8293068783069,28.6731164021164,29.5169259259259,30.3607354497354,31.2045449735450,32.0483544973545,32.8921640211640,33.7359735449736,34.5797830687831,35.4235925925926,36.2674021164021,37.1112116402116,37.9550211640212,38.7988306878307,39.6426402116402,40.4864497354497,41.3302592592593,42.1740687830688,43.0178783068783,43.8616878306878,44.7054973544974,45.5493068783069,46.3931164021164,47.2369259259259,48.0807354497355,48.9245449735450,49.7683544973545,50.6121640211640,51.4559735449735,52.2997830687831,53.1435925925926,53.9874021164021];
else % user input no TIRS instrument choice or number <> 1, 2, coincides with original estimates
    TIRS.fl = 8020; % focal length of optics in microns, spatial dimension
    TIRS.fl2 = 6471; % focal length of optics in microns, spectral dimension
    TIRS.R = 3600*ones(TIRS.spectral_channels,TIRS.spatial_scenes); % latest CBE 2/25/2021
    TIRS.R(1:4,:) = (4/16)*TIRS.R(1:4,:); % 4 thermopiles instead of 16, in rows zero through four for zero-order
    TIRS.Dstar = TIRS.R*1e-4*TIRS.d*sqrt(TIRS.tau)/TIRS.noise;
    TIRS.x_fp = 0.0; %model values
    TIRS.tilt_fp = deg2rad(0.0); %tilt of focal plane, rotation about x-axis
    TIRS.y_fp = 0.0;%50.0; % offset of focal plane in spatial dimension
    TIRS.arb = ones([64,8]);
    TIRS.centers = [0,0.843809523809524,1.68761904761905,2.53142857142857,3.37523809523810,4.21904761904762,5.06285714285714,5.90666666666667,6.75047619047619,7.59428571428571,8.43809523809524,9.28190476190476,10.1257142857143,10.9695238095238,11.8133333333333,12.6571428571429,13.5009523809524,14.3447619047619,15.1885714285714,16.0323809523810,16.8761904761905,17.7200000000000,18.5638095238095,19.4076190476190,20.2514285714286,21.0952380952381,21.9390476190476,22.7828571428571,23.6266666666667,24.4704761904762,25.3142857142857,26.1580952380952,27.0019047619048,27.8457142857143,28.6895238095238,29.5333333333333,30.3771428571429,31.2209523809524,32.0647619047619,32.9085714285714,33.7523809523810,34.5961904761905,35.4400000000000,36.2838095238095,37.1276190476191,37.9714285714286,38.8152380952381,39.6590476190476,40.5028571428571,41.3466666666667,42.1904761904762,43.0342857142857,43.8780952380952,44.7219047619048,45.5657142857143,46.4095238095238,47.2533333333333,48.0971428571429,48.9409523809524,49.7847619047619,50.6285714285714,51.4723809523810,52.3161904761905,53.1600000000000];
end

% values are 5000-6000 V/W for regular pixels 2000-3000 V/W in band zero
% these may be optimistic because I am assuming noise is only 120 nV/Hz^1/2
% however, the source spreadsheet gives noise in counts/Hz^1/2 and values
% are near unity, which is consistent with 1/TIRS.g ~ 131 nV/Hz^1/2

%% definition of focal plane offsets, with respect to coordinate system origin in center of ideal focal plane

TIRS.z_fp = 0.0; % offset of physical focal plane coordinate
TIRS.tau_fp = deg2rad(0.0);%0.25; %clocking error of focal plane, rotation about z-axis
TIRS.tip_fp = deg2rad(0.0);%2.9; %tip of focal plane, rotation about y-axis

% these are the spectrometer design parameters and derived quantities
% f#'s to 4 significant figures provide by B. Cameron memo #9
TIRS.fno = 2.005; % f# of the optics in the spatial
TIRS.fno2 = 1.618; % spectral f# 
% focal lengths to 3 significant figures provide by B. Cameron memo #9
TIRS.fl_nom = 8020; %nominal focal length for scaling image using measured focal length, spatial dimension
TIRS.fl2_nom = 6471; %nominal focal length for scaling image using measured focal length, spectral dimension

TIRS.bandwidth = 53.16; %microns
TIRS.dispersion = TIRS.d*(TIRS.spectral_channels-1)/TIRS.bandwidth;%210;%ideal dispersion spectral microns/spatial microns
TIRS.sw = slit_flag*TIRS.d; %flag is ratio of slit width to pixel size
TIRS.channel_sampling = TIRS.bandwidth/(TIRS.spectral_channels-1);% the first channel is non-spectral, or zero order
TIRS.channel_center_wavelens = (0:TIRS.spectral_channels-1)'*TIRS.channel_sampling; %now zero-based!
TIRS.sl = 3940;% slit length in microns

%% definition of slit offsets, with respect to ideal slit position
%% asuuming slit is perfect rectangle with non-ideal placement/clocking 
TIRS.x_slit = 0.0; % offset of slit in spectral dimension
TIRS.y_slit = 0.0; % offset of slit in spatial dimension
TIRS.tau_slit = 0.0;% deg2rad(1.0); % clocking error of slit B. testing tolerance 7/21/2022 this breaks really good! 9/26/22
TIRS.slit_image_position = -5670; %offset between ideal relative positions of slit and focal plane in microns

%% 0th order smile, to be applied geometrically at projection to surface
%% degrees, scenes 1-8, @-13.93 -9.90 -5.92 -1.97 1.97 5.29 9.90 13.93
TIRS.smile = deg2rad([1.265 0.627 0.219 0.021 0.021 0.219 0.627 1.265]);
%set up optical (view) factors that account for thermal zones, same
%definitions as in BRF model
%uncertain factors here
TIRS.theta = 2*atan(0.5/TIRS.fno); %factor of two???, this is the cone angle 
TIRS.theta2 =  2*atan(0.5/TIRS.fno2);
TIRS.A_inFOV=pi*4*(sin(TIRS.theta/2))^2; %solid angle in steradians
TIRS.A_inFOV=pi*TIRS.A_inFOV*4; % extra factor of four for unknown reason - allows W/m2 result to convert A_all back to temperature input with stefan boltzmann constant
TIRS.A_all=pi*4;%(sin(2*atan(0.5/0.0001)))^2; %solid angle in steradians
TIRS.A_all = TIRS.A_all*4; % extra factor of four for unknown reason - allows W/m2 result to convert back to temperature input with stefan boltzmann constant

TIRS.reflectance_eff = 0.75;
%potential for spectral dependent reflective model (not yet implemented)
Aluminum.Reflectivity = [.926 .920 .916 .907 .868 .890 .940 .980 .987 .990 .994];
Aluminum.Wavelengths =  [.248 .400 .532 .633 .800 .900 1.00 3.00 10.6 20.0 100.];
TIRS.Reflectivity = interp1(Aluminum.Wavelengths,Aluminum.Reflectivity,TIRS.channel_center_wavelens);
%first channel (zero order) give NaN, so replace with average
TIRS.Reflectivity(1) = sum(TIRS.Reflectivity(2:TIRS.spectral_channels))/(TIRS.spectral_channels-1);
TIRS.num_optics = 6;
TIRS.Reflectivity = TIRS.Reflectivity.^TIRS.num_optics;

% filter locations
TIRS.filter_channel = ...
    [  1  1  0  0  2  2  2  2  0  0  3  3  3  3  3  3  3  0  0  4];
%      0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19
TIRS.filter_channel = [TIRS.filter_channel, ...
       4  4  4  4  4  4  4  4  4  4  4  4  4  4  4  0  0  5  5  5];
%     20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39
TIRS.filter_channel = [TIRS.filter_channel, ...
       5  5  5  5  5  5  5  5  5  5  5  5  5  5  5  5  5  5  5  5];
%     40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59
TIRS.filter_channel = [TIRS.filter_channel, 5, 5, 5, 5];
%                                           60 61 62 63

end
