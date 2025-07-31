% remember to add path for the parent dir, for the exitance
% function. e.g.:
%
% addpath /home/mmm/projects/TIRS_spectral_model

% loads SRF, shaped [nwave, nchannel], and 
% wavelen, shaped [nwave].

[script_path, ~, ~] = fileparts(which('PREFIRE_NEDT_model_v5'));
orig_path = path;
path(orig_path, fullfile(script_path, 'functions'));

pd = configure_toplevel_IO;  % Get/set some shared filepaths and directories

%fpath_to_load = fullfile(pd.ancillary_data_dir, 'SRF', ...
%                             'PREFIRE_SRF_v0.11_360_2022-11-03.mat'); % inst 1
fpath_to_load = fullfile(pd.ancillary_data_dir, 'SRF', ...
                 'PREFIRE_SRF_v0.11_360_2022-11-15_instrument2.mat'); % inst 2
load(fpath_to_load);


% in March 28 2021 version, the channel wavelengths became scene dependent
%% Instrument parameters
TIRS = prefire_TIRS(2,1,1); %change the third input for TIRS choice, 0 = model, 1 = FPA101, 2 = FPA102

%version 5 role in empirical correction from radcal
for i=1:length(wavelen)
    SRF(i,:,:)= reshape(SRF(i,:,:),[64,8]).*TIRS.arb;
end

TIRS.NEP = TIRS.noise./TIRS.R;
%get Dstar from prefire_TIRS now
%TIRS.Dstar = sqrt(1e-4*TIRS.d*1e-4*TIRS.d*TIRS.tau)./TIRS.NEP; %convert to cm/W 


%TIRS.Dstar(1:TIRS.spectral_channels,1:TIRS.spatial_scenes) = 0.5e9; %minimum acceptable value (changed from 0.7 to 0.5 on 2/4/2021)

TIRS.NEP = sqrt(1e-4*TIRS.d*1e-4*TIRS.d*TIRS.tau)./TIRS.Dstar; %recalculate from minimum acceptable value

TIRS.NEE=(TIRS.NEP*4*TIRS.fno*TIRS.fno)/(1e-4*TIRS.d*1e-4*TIRS.d); %Noise Equivalent Exitance from Scene in W/cm2 
%fix for anamorphic FOVs
TIRS.NEE=(TIRS.NEP*4*TIRS.fno*TIRS.fno2)/(1e-4*TIRS.d*1e-4*TIRS.d); %Noise Equivalent Exitance from Scene in W/cm2 
% reduce for constricted FOV, linear correction 
TIRS.NEE = TIRS.NEE*TIRS.theta2/TIRS.theta;


TT=50:50:500; %target temperature in K  TT=3:0.1:500;

[num_wave, num_channel] = size(SRF);

%% start loop over T, etc.
T = (200:10:300)';
num_T = size(T,1);
SNR = zeros([TIRS.spectral_channels, num_T, TIRS.spatial_scenes]);
dw = wavelen(2) - wavelen(1);
ex = zeros([num_wave, num_T]);

% compute blackbody spectrum on same grid as SRF.
for n = 1:num_wave;
    start = wavelen(n) - 0.5 * dw;
    %start = abs(wavelen(n) - 0.05);
    stop = wavelen(n) + 0.5 * dw;
    %stop = wavelen(n) + 0.05;
    for t = 1:num_T;
        ex(n,t) = exitance(start, stop, T(t));
        %ex(n,t) =  10*exitance(start,stop,T(t))/pi; %same expression and validated in Emission PreCalculator, units of cm-2 match NEE?
    end
end

% to compute channel SNR, just take dot product with the blackbody
% spectrum. Note the ex array already represents an integrated
% quantity so we don't need to multiply by dw (I think)

for k = 1:TIRS.spectral_channels
    for j = 1:TIRS.spatial_scenes
        for t = 1:num_T
            SNR(k,t,j) = (ex(:,t)' * SRF(:,k,j)) ./ TIRS.NEE(k,j);    
        end
    end
end

%computing/plotting just first spatial scene
TD = T(7) - T(5);
for j = 1:TIRS.spatial_scenes
    NETD(:,j) = TD./(SNR(:,7,j)-SNR(:,5,j));
end

channel_wave_centers = mean(channel_wavelens,2);


figure()
plot(channel_wave_centers, reshape(SNR(:,6,:),[64,8]));
xlabel('Wavelength [um]')
ylabel('SNR')

figure()
plot(channel_wave_centers, NETD);
xlabel('Wavelength [um]')
ylabel('NETD')
%save("instrument_1_NETD_NEDR_SRF.mat",'NEDR','NETD','SRF')
save("instrument_2_NETD_NEDR_SRF.mat",'NEDR','NETD','SRF')

path(orig_path);
