function [geometry] = geometry_read2(varargin)% 'filename'
%% PREFIRE geometry data Reader 
% 04062020 B.
if(length(varargin)>0)
    filename=varargin{1};
else
    [filename, pathname]=uigetfile('*.bin','select geometry file'); %default is to prompt for file
end

%% Initialize the Workspace and Preallocate the Data
%close all %close any open windows
%clear   %initialize the workspace

fdata=dir(filename);

%% Read the File
fi=fopen(filename,'rb');
geometry.slit_width=fread(fi,1,'int16');%'uint16');
geometry.pix_width=fread(fi,1,'int16');%'uint16');
geometry.num_spectral=fread(fi,1,'int16');%'uint16');
geometry.num_spatial=fread(fi,1,'int16');%'uint16');
geometry.x_slit=fread(fi,1,'int16');%'uint16');%-32768;
geometry.y_slit=fread(fi,1,'int16');%'uint16');%-32768;
geometry.tau_slit=fread(fi,1,'int16');%'uint16');%-32768;
geometry.x_fp=fread(fi,1,'int16');%'uint16');%-32768;
geometry.y_fp=fread(fi,1,'int16');%'uint16');%-32768;
geometry.z_fp=fread(fi,1,'int16');%'uint16');%-32768;
geometry.tau_fp=fread(fi,1,'int16');%'uint16');%-32768;
geometry.tip_fp=fread(fi,1,'int16');%'uint16');%-32768;
geometry.tilt_fp=fread(fi,1,'int16');%'uint16');%-32768;
geometry.xy_p = zeros([2,geometry.num_spectral,geometry.num_spatial,5]);
%xy_p = fread(fi,[2,num_spectral,num_spatial,5],'uint16')-32768;
for m = 1:5
    for j = 1:geometry.num_spatial
        for k = 1:geometry.num_spectral
            geometry.xy_p(1,k,j,m) = fread(fi,1,'int16');%'uint16');%-32768; 
            geometry.xy_p(2,k,j,m) = fread(fi,1,'int16');%'uint16');%-32768; 
        end
    end
end
geometry.slit=fread(fi,[2,5],'int16');%'uint16');%-32768;

fclose(fi);

end
