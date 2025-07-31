function [L1b, filepath] = prefire_read_full_spectrum(varargin)
%% PREFIRE L1b reader modified for full spectrum PCRTM output
% update based on newer sims - this is now 'Anc-FullSpectrum'
% Ancillary product, with Geometry and full spectra in the same
% product file. (As of May 2022.)

%% Initialize the Workspace and Preallocate the Data
close all %close any open windows
%clear   %initialize the workspace

if(length(varargin)>0)
    plot_flag=varargin{1};
else
    plot_flag=0; %default value
end

if(length(varargin)>1)
    radfile = varargin{2};
    filepath = radfile;
else
    [radfile, pathname]=uigetfile('*.nc','select PREFIRE_SAT1_1B-RAD file'); %default is to prompt for file
    filepath = append(pathname,radfile);
end

L1b.latitudes = ncread(filepath, '/Geometry/latitude');
L1b.longitudes = ncread(filepath, '/Geometry/longitude');
L1b.times = ncread(filepath, '/Geometry/time_UTC');
L1b.seconds(:) = ((L1b.times(3,:)*24+L1b.times(4,:))*60+L1b.times(5,:))*60+L1b.times(6,:)+L1b.times(7,:)/1000; % seconds

L1b.full_spec_rad = ncread(filepath, '/Anc-Fullspectra/radiance');
%will need to pull in the wavelength grid too 
L1b.wavenum = ncread(filepath, '/Anc-Fullspectra/wavenum');


if plot_flag
    hires_wavelens = 10000./L1b.wavenum;
    test = 10.*L1b.full_spec_rad(:,1,1)./hires_wavelens./hires_wavelens;
    figure
    plot(hires_wavelens,test,wavelens,ans.spec_rad(:,1,1))
%    plot_L1b(filepath,L1b.latitudes,L1b.longitudes,L1b.seconds,L1b.spec_rad);
end

end
