function [housekeeping] = get_housekeeping(latitudes);
%this function determines the internal temperatures of TIRS based on orbit position

% old definitions
% phases(1) = deg2rad(9.0); %Filter 
% phases(2) = deg2rad(10.0); %Detector
% phases(3) = deg2rad(0.0); %Spacecraft
% phases(4) = deg2rad(3.0); %Calibration target
% phases(5) = deg2rad(9.1); %Filter'
% phases(6) = deg2rad(8.0); %Grating
% phases(7) = deg2rad(7.0); %Toroid
% phases(8) = deg2rad(10.1); %Detector'

% base_temperature = 280.0; %Kelvin (280 K cold case, 300 K hot case)
% amplitudes(1) = 3.0; %Kelvin (smaller in hot case) Filter
% amplitudes(2) = 3.0; %Kelvin (smaller in hot case) Detector
% amplitudes(3) = 10.0; %Kelvin (smaller in hot case) Spacecraft
% amplitudes(4) = 3.0; %Kelvin (smaller in hot case) Calibration Target
% amplitudes(5) = 3.0; %Kelvin (smaller in hot case) Filter'
% amplitudes(6) = 5.0; %Kelvin (smaller in hot case) Grating
% amplitudes(7) = 5.0; %Kelvin (smaller in hot case) Toroid
% amplitudes(8) = 3.0; %Kelvin (smaller in hot case) Detector'

% new definitions
phases(1) =  7.0; % Toroidal mirror
phases(2) =  0.0; % telescope - assume it maps to old 'spacecraft'
phases(3) =  9.0; % filter
phases(4) =  3.0; % Calibration target
phases(5) =  3.0; % Calibration target
phases(6) =  3.0; % Calibration target
phases(7) = 10.0; % FPA center/back (detector?)
phases(8) =  9.6; % FPA mounting bracket - assume midway between FPA
                  % center and filter?
phases = deg2rad(phases);

base_temperature = 280.0;
amplitudes(1) =  5.0;
amplitudes(2) = 10.0;
amplitudes(3) =  3.0;
amplitudes(4) =  3.0;
amplitudes(5) =  3.0;
amplitudes(6) =  3.0;
amplitudes(7) =  3.0;
amplitudes(8) =  3.0;

housekeeping = zeros([8,size(latitudes,2)]);
for m=1:size(latitudes,2)
    housekeeping(:,m) = base_temperature + ...
        amplitudes.*(1+cos(deg2rad(latitudes(1,m))+phases));;
end

end
