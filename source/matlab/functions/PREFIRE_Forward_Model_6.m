function PREFIRE_Forward_Model_6(filepath,varargin)
%   ncdisp("C:/Users/yyy/Documents/My Work/PREFIRE/PREFIRE_L1/03UTC/PREFIRE_SAT1_1B-RAD_S00_R00_20210419000000_00001.nc")
% This is the sixth iteration of B. forward model for PREFIRE
% the sixth version is the first with hi-res data from PCRTM available in
% orbit sampled aggregates
%
% required inputs:
% filepath - "Anc-FullSpectra" netCDF product file, containing
%   hyperspectral simulated IR data and sensor Geometry group.
%
% optional inputs:
% plot_flag = 0 or 1, depending on whether you want the diagnostic
%   plots created (default 0).
% rng_seed = 'none', or some integer.
%   if set to 'none', then the rng is not initialized and the rng
%   will proceed to make a new random number stream for the
%   simulated detector noise.
%   set this to some specific integer to make identical noise
%   realizations for testing. (default 'none')
% add_noise = 0 or 1, depending on whether you want to add detector
%   noise to the simulated counts. This was added for testing, to
%   see the effect of the detector quantization. Normally you want
%   to leave this on (default 1).

pd = configure_toplevel_IO;  % Get/set some shared filepaths and directories

if length(varargin) >= 3
    add_noise = varargin{3};
else
    add_noise = 1;
end

if length(varargin) >= 2
    rng_seed = varargin{2};
else
    rng_seed = 'none';
end

if length(varargin) >= 1
    plot_flag = varargin{1};
else
    plot_flag = 0;
end

% Start by importing in a simulated L1b file
% plot flag and path\filename may be passed, if blank then user is prompted for filename
[L1b, filepath] = prefire_read_full_spectrum(0,filepath);
%now have L1b.full_spec_rad and L1b.wavenum that are 'parents' to L1b.spec_rad
%L1b.spec_rad = fillmissing(L1b.spec_rad,'constant',0); % there were some NaNs in at least 1 PCRTM output
% This routine will hold all design parameters for instrument
[TIRS] = prefire_TIRS(2,1); %first argument is slit/pixel ratio, second is low/high gain (zero/non-zero)

% Set up instrument thermal data, need only one latitude 'strip' choose scene 4,
% this is the instrument thermal model
housekeeping = get_housekeeping(L1b.latitudes(4,:)); %may want to add beta angle at some point
fpath_to_load = fullfile(pd.ancillary_data_dir, 'BRF', ...
                                           'PREFIRE_BRF_v0.05_2020-05-31.mat');
load(fpath_to_load); %file should be matched to slit width and TIRS #
T_zero = T_cal(1)-1; % origin defined in PREFIRE Emission Precalculator

% Set up calibration cadence
cal_indices = get_cal_indices(L1b.seconds,TIRS.cadence);

% prepare local variables
num_images = [TIRS.spatial_scenes, size(L1b.latitudes,2)];
sig = zeros(TIRS.spectral_channels,num_images(1),num_images(2));

% versions < 6 fill signal vector with L1b data and its sum
%sig(1,:,:)= reshape(sum(L1b.spec_rad,1),num_images);
%sig(2:TIRS.spectral_channels,:,:) = L1b.spec_rad;%*A_inFOV2; 

%version 6, convolute L1b.full_spec_rad with instrument function

% matfile contains 'wavelen', size [n_wave, 1] and 'SRF', shape
% [n_wave, n_channel, n_scene].
fpath_to_load = fullfile(pd.ancillary_data_dir, 'SRF', ...
                                     'PREFIRE_SRF_v0.10.4_360_2021-03-27.mat');
load(fpath_to_load); %B. update 1/24/22 to slightly newer version (3-15 ->0 3-27)

% compute per-wavenum conversion factor to convert the simulated
% full spectra to per-wavelength.
% conversion includes mW->W, and cm->um, resulting in factor of 10.
L1b_wavelen = 1e4 ./ L1b.wavenum;
wave_conversion = 10.0 ./ L1b_wavelen.^2;
% using broadcasting here. L1b.full_spec_rad is (n_wave, n_xtrack, n_atrack)
full_spec_wavelen = L1b.full_spec_rad .* wave_conversion;

% now, interpolate the full_spec to the wavelen grid for SRF.
% this seems to do the right thing, broadcasting over axis 2 and 3
full_spec_wavelen_interp = interp1(L1b_wavelen, full_spec_wavelen, wavelen, 'linear', 0.0);

%ideally this could be completely vectorized
%sig = sum(SRFinterp.*L1b.full_spec_rad,1); %64 sums over N(hires), get N(images)x64x8
%for now this channel by channel, scene by scene appears to be working
SRFsums = squeeze(sum(SRF, 1));
for k=1:64
    for j=1:8
        if SRFsums(k,j) > 0
            localSRF = SRF(:,k,j);
            localSpec = full_spec_wavelen_interp(:,j,:);
            sig(k,j,:) = reshape(sum(localSRF.*localSpec,1),[1, length(L1b.seconds)]) / SRFsums(k,j);
        end
    end % for j
end %for k

%B. add zero order sum 11/22/21
       %overwrite zero order values
       %zeroth order space and calibration values are sums of all channels
       %(not quite true with SRF gaps at blocking between filters)
       %there is no grating efficiency or reflection coefficient applied!
for j=1:8
    sig(1,j,:) = sum(sig(2:64,j,:));
end

% interlace signal vector with calibration data
%encoder=zeros(num_intrack,1);
encoder = ones(size(L1b.latitudes,2),1)*TIRS.nadir_enc; % nadir pointing encoder value
for mm=1:size(cal_indices,2)
   if (cal_indices(mm)+3 < size(L1b.latitudes,2)) % prevents writing past end of array if cal cadence places a sequence within last few images of sim data
       T_cal = round(housekeeping(4,cal_indices(mm)))-T_zero; % get calibration target temperature and shift to precalculated grid counter
       encoder(cal_indices(mm)) = TIRS.space_enc;
       encoder(cal_indices(mm)+1) = TIRS.space_enc;
       encoder(cal_indices(mm)+2) = TIRS.cal_enc;
       encoder(cal_indices(mm)+3) = TIRS.cal_enc;
       %spectral channels
       sig(:,:,cal_indices(mm)) = SRFMS(:,:);% space scenes convolved with SRFs, no zero order
       sig(:,:,cal_indices(mm)+1) = SRFMS(:,:);% space scenes convolved with SRFs, no zero order
       sig(:,:,cal_indices(mm)+2) = SRFMC(:,:,T_cal);% calibration target scenes convolved with SRFs, no zero order
       sig(:,:,cal_indices(mm)+3) = SRFMC(:,:,T_cal);% calibration target scenes convolved with SRFs, no zero order
       %overwrite zero order values
       %zeroth order space and calibration values are sums of all channels
       %(not quite true with SRF gaps at blocking between filters)
       %there is no grating efficiency or reflection coefficient applied!
       sig(1,:,cal_indices(mm)) = sum(SRFMS(:,:),1);
       sig(1,:,cal_indices(mm)+1) = sum(SRFMS(:,:),1);
       sig(1,:,cal_indices(mm)+2) = sum(SRFMC(:,:,T_cal),1);
       sig(1,:,cal_indices(mm)+3) = sum(SRFMC(:,:,T_cal),1);
   end %if cal_indices
end    
%this model has no geometric signal scaling factors, it currently assumes
%all of scene contributes to thermopile signal

% background signal is everything that is not scene and is difference
% between that scene 'temperature' and detector 'temperature'
scene = TIRS.A_inFOV/TIRS.A_all;
not_scene = 1-scene;

background = BRFMD(:,round(housekeeping(2,:))-T_zero); % filter reflection minus detector emission at detector temperature
background = background + BRFMF(:,round(housekeeping(1,:))-T_zero); %add on filter emissions at filter temperature
background = background + BRFMD(:,round(housekeeping(7,:))-T_zero); %add on instrument emissions at toroid temperature

TIRS_length = size(L1b.latitudes,2);

for k=1:64
    f = TIRS.filter_channel(k);
    if (f > 0)
        local_background = not_scene*background(f,:);
    else
        local_background = zeros([1,TIRS_length]); %these are the blocking filters
    end
    for j=1:8
        temp = scene*sig(k,j,:);
        sig(k,j,:) = reshape(temp,[1,TIRS_length])+local_background;
    end
end

if add_noise
    % create random noise vector for signal
    if rng_seed ~= 'none'
        rng(rng_seed)
    end
    pd = makedist('Normal');
    r = random(pd,[TIRS.spectral_channels,num_images(1),num_images(2)]); %looks like this has a std deviation of 2 distributed around 0
    r = r*TIRS.noise; % 60 nV for high gain, 90 nV for low gain
else
    r = zeros([TIRS.spectral_channels,num_images(1),num_images(2)]);
end

d = TIRS.d/1e6; % convert from microns to meters
for m=1:num_images(2)
    sig(:,:,m) = TIRS.R.*sig(:,:,m)*d^2;
end

if plot_flag
    test = reshape(sig(1,1,:),[1,TIRS_length]);
    test2 = reshape(sig(11,1,:),[1,TIRS_length]);
    test3 = reshape(sig(21,1,:),[1,TIRS_length]);
    test4 = reshape(sig(31,1,:),[1,TIRS_length]);
    test5 = reshape(sig(41,1,:),[1,TIRS_length]);
    test6 = reshape(sig(51,1,:),[1,TIRS_length]);
    plot(1:TIRS_length,scene*test,1:TIRS_length,scene*test2,1:TIRS_length,scene*test3,1:TIRS_length,scene*test4,1:TIRS_length,scene*test5);
end
    
sig = TIRS.o+TIRS.g*(r+sig);

% create a rounded off counts value here before plotting!
sig = round(sig,0);
if plot_flag 
    forward_model_plots(TIRS.spectral_channels, num_images, TIRS.channel_center_wavelens, L1b.latitudes, sig, '');
end    


signew = ROIC_map_oi(size(L1b.latitudes),TIRS.spectral_channels,sig); % re-order data to assumed readout format
% output data3 now carries cadence info instead of determining cadence from
% num, which is size(cal_indices,2)-1 here, it is only used for determining
% a sequence number
% TIRS.spectral_channels may be one index larger than output_data2 is expecting

%From ST-ICD:
%Level-0 files are named with the following convention:

%prefire_01_payload_<yyyy>_<mm>_<dd>_<yyyy>_<mm>_<dd>.csv
%prefire_02_payload_<yyyy>_<mm>_<dd>_<yyyy>_<mm>_<dd>.csv
[input_dir, input_file, input_ext] = fileparts(filepath);
fileout = ['PREFIRE_L0_data/', input_file, '.bin'];
fileout = strrep(fileout,'_Anc-Fullspectra','_0-RAD');
disp(fileout);
output_data3(fileout,num_images,TIRS.cadence,L1b.seconds,TIRS.spectral_channels,housekeeping,encoder,signew);
% use inner two footprints as a geodetic ground point and create ECI s/c
% coordinates by setting altitude, then write all of the fields in the
% COSMOS output
fileout = ['PREFIRE_L0_data/', input_file, '.csv'];
fileout = strrep(fileout,'_Anc-Fullspectra','_0-GPS');
disp(fileout);
csvID = fopen(fileout,'w');
string = sprintf('%s\n', ...
    ['time,ATTITUDE_VALID_val,BATTERY1_TEMP_val,BATTERY2_TEMP_val,BATTERY_CURRENT_val,BATTERY_VOLTAGE_val,'...
     'BOX1_TEMP_val,BUS_VOLTAGE_val,CADET_RECEIVER_TEMP_val,CADET_TEMP_val,GPS_VALID_val,HEATER_STATUS1_val,'...
     'HEATER_STATUS2_val,HEATER_STATUS3_val,IO11_IN_SEP_MON_val,IO12_IN_BAT2_CHG_STAT_val,'...
     'IO16_IN_BAT1_CHG_STAT_val,IO17_OUT_5V_val,IO19_OUT_RELEASE2_val,IO1_OUT_CADET_PA_val,'...
     'IO20_OUT_REL_REG_val,IO21_OUT_PL_5V_val,IO22_OUT_PL_7V5_val,IO23_OUT_7V5_val,IO24_OUT_RELEASE1_val,'...
     'IO25_OUT_OUT_RELEASE3_val,IO26_OUT_PL_15V_val,IO2_OUT_GSTAR_val,IO3_OUT_GPS_val,IO4_OUT_PL_HTR_val,'...
     'IO5_OUT_CADET_val,IO6_OUT_BUS_HTR_val,IO7_OUT_BUS_HTR_val,IO8_OUT_SM_HTR_val,MEAS_WHEEL_CURRENT1_val,'...
     'MEAS_WHEEL_CURRENT2_val,MEAS_WHEEL_CURRENT3_val,MOTOR1_TEMP_val,MOTOR2_TEMP_val,MOTOR3_TEMP_val,'...
     'POSITION_ECEF1_val,POSITION_ECEF2_val,POSITION_ECEF3_val,POSITION_WRT_ECI1_val,POSITION_WRT_ECI2_val,'...
     'POSITION_WRT_ECI3_val,Q_BODY_WRT_ECI1_val,Q_BODY_WRT_ECI2_val,Q_BODY_WRT_ECI3_val,Q_BODY_WRT_ECI4_val,'...
     'TAI_SECONDS_val,TIME_VALID_val,TORQUE_ROD_ENABLE1_val,TORQUE_ROD_ENABLE2_val,TORQUE_ROD_ENABLE3_val,'...
     'TRACKER2_OPERATING_MODE_val,TRACKER_OPERATING_MODE_val,VELOCITY_WRT_ECI1_val,VELOCITY_WRT_ECI2_val,'...
     'VELOCITY_WRT_ECI3_val,VOLTAGE_12P0_val,VOLTAGE_3P3_val,VOLTAGE_5P0_val,VOLTAGE_8P0_val']);
fwrite(csvID,string);
for i=1:20:TIRS_length
    if i < (TIRS_length-1)
        iprime = i + 1;
    else
        iprime = i;
        i = i - 1;
    end
    lat = (L1b.latitudes(4,i)+L1b.latitudes(5,i))/2;
    lon = (L1b.longitudes(4,i)+L1b.longitudes(5,i))/2;
    alt = 525000;
    ECI = lla2eci([lat lon alt],reshape(L1b.times(1:6,i),[1,6]))/1000;
    latp = (L1b.latitudes(4,iprime)+L1b.latitudes(5,iprime))/2;
    lonp = (L1b.longitudes(4,iprime)+L1b.longitudes(5,iprime))/2;
    altp = 525000;
    ECIp = lla2eci([latp lonp altp],reshape(L1b.times(1:6,iprime),[1,6]))/1000;
    Qt = norm(ECI-ECIp)/TIRS.tau; % velocity in km/s
    Q = (ECIp-ECI)./Qt; %unit velocity vector
    Qt = Qt/2;
    string = sprintf( ...
        ['"%d-%d-%dT%02d:%02d:%02d.%03dZ","YES","","","","","","","","","NO","OFF","OFF","OFF",' ...
         '"","","","","","","","","","","","","","","","","","","","","","","","","","","","","",' ...
         '%f,%f,%f,%f,%f,%f,%f,%f,'...
         '"","","","","","","","","","","","",""\n'], ...
	 L1b.times(1,i),L1b.times(2,i),L1b.times(3,i),L1b.times(4,i), ...
	 L1b.times(5,i),L1b.times(6,i),L1b.times(7,i),ECI(1),ECI(2),ECI(3), ...
	 Qt,Q(1),Q(2),Q(3),L1b.seconds(1,i));
    fwrite(csvID,string);
end   
fclose(csvID);

end
