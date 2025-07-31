function [gm, wl] = geometry_v3(varargin)
% geometric model of PREFIRE TIRS
%
% input arguments:
%   TIRS number: integer equal to 1 or 2, specifying which instrument to
%       simulate. required.
%   doplots: logical to specify if plots should be made.
%       this is optional, defaults to false if inspecified.
%
% 10/17/22, still having trouble with plot routines, rewrite i/o again,
% this time to allow matlab i/o
% minor updates 9/26/22 for reading back in goemetry.bin and wavelength.bin
% for passing structures to plot routines
%
% on 2/22/2021 it became apparent that the anamorphic nature of the optics
% will affect the defocus calculated in this code, therefore the two
% dimensions are now specifically addressed in geometry_v3.m

% on 5/29/20 B. created prefire_TIRS.m to hold much of the input values,
% geometry_v2 has been modified to call that one...
% B. 3/19/2020-3/24/2020
% change values for x_fp, y_fp, z_fp, tau_fp, tip_fp, tilt_fp, x_slit,
% y_slit and tau_slit to obtain different geometric and spectrographic
% instrument results, 
% geometric output in xy_p arrays has effective image for each
% pixel in microns, which projects onto ground pixels with known orbital
% height, throughput is not considered, so darkening of pixels not
% illuminated by the slit is TBD, but geometric displacements and 
% distortions should be accurate.
% xy_s array is currently unused 
% the wavelength arrays compare center and edge pixel positions vs. ideal

TIRS_number = varargin{1};
if nargin > 1
    doplots = varargin{2};
else
    doplots = false;
end

pd = configure_toplevel_IO;

% This routine will hold all design parameters for instrument
[TIRS] = prefire_TIRS(0.7007, pd, 2, 0, TIRS_number); %3rd argument is slit/pixel ratio, 4th is low/high gain (non-zero/zero)

%% instrument definitions
d = TIRS.d; % pixel pitch in microns
d2 = TIRS.d2; % pixel strip gap in microns
slit_width = TIRS.sw;% slit width in microns
slit_length = TIRS.sl;% slit length in microns
dispersion = TIRS.dispersion;%ideal dispersion spectral microns/spatial microns

%% definition of array format
num_spectral = TIRS.spectral_channels; % pixel count in x, or spectral, coordinate
num_spatial = TIRS.spatial_scenes; % pixel count in y, or spatial, coordinat

%% optical parameters
theta_spatial = atan(0.5/TIRS.fno); %beam divergence angle, measured from normal about 14 degrees
theta_spectral = atan(0.5/TIRS.fno2); % allows for anamorphic defocus

%% Filter defocus and focus spread information from 
%% PREFIRE-OXF-SP-002 2-0 PREFIRE substrate thicknesses and defocus.xlsx
%% received from Simon Calcutt on 4/12/2020
% not implemented in algorithm yet
%Filter |   Total Focus Offset (microns) | Focus offset variance (microns) [spread]
%_________________________________________________________________
%0      |      0                          |   25
%1      |    -57                          |   32
%2      |      0                          |   63
%3      |      -                          |   16
%4      |    -37                          |   60
% ZnS (F#1) and CdTe (F#3) are 'real' whereas other filters are TBC

%% Current requirements for detector positioning wrt true optical focal
%% plane per Andrew Houck communication 4/29/2020
% max_tip = 17 mrad (spatial direction) 0.974 degrees (enter into tilt_fp)
% max_tilt = 50 mrad (spectral direction) 2.865 degrees (enter into tip_fp)
%% does not include as-built parallelism of true focal plane wrt to grating plane
% in-plane requirements 
% max_x = 0.05 mm (decenter)
% may_y = 0.05 mm (decenter)
% R_z = 4.2 mrad (clocking) 0.241 degrees (enter into tau_fp)
%% B. notes that Houck definition of tip/tilt seems to be opposite of B. 
%% definition, so these are reversed upon entry below
%% B. also notes that no translation in z (out of focus) was provided

%% definition of slit offsets, with respect to ideal slit position
%% asuuming slit is perfect rectangle with non-ideal placement/clocking 
x_slit = TIRS.x_slit; % offset of slit in spectral dimension
y_slit = TIRS.y_slit; % offset of slit in spatial dimension
tau_slit = TIRS.tau_slit; % clocking error of slit

%% definition of focal plane offsets, with respect to coordinate system origin in center of ideal focal plane
slit_image_position = TIRS.slit_image_position; %offset between ideal relative positions of slit and focal plane in microns
x_fp = TIRS.x_fp;%50.0; % offset of focal plane in spectral dimension
y_fp = TIRS.y_fp;%50.0; % offset of focal plane in spatial dimension
z_fp = TIRS.z_fp; % offset of physical focal plane coordinate
tau_fp = TIRS.tau_fp;%0.25; %clocking error of focal plane, rotation about z-axis
tip_fp = TIRS.tip_fp;%2.9; %tip of focal plane, rotation about y-axis
tilt_fp = TIRS.tilt_fp;%1.0; %tilt of focal plane, rotation about x-axis

%% corner definitions
xyc = zeros([2,4]);
for m=1:4
    xyc(1,m) = 1 - 2*mod(m,2); % alternates +,+,-,-      
    xyc(2,m) = 1 - 2*floor((m-1)/2); % alternates +,-,+,-
    
end

%% definition of image positions within slit, ideal initialization
xy_s = zeros*[2,num_spatial,5]; %512 x,y coordinates for center, four corners
slit = zeros([2,5]);
slit(1,1) = 0.0; %translate to slit_image_position later
slit(2,1) = 0.0;
for m=1:4
   slit(1,m+1)=slit(1,1)+xyc(1,m)*slit_width/2; 
   slit(2,m+1)=slit(2,1)+xyc(2,m)*slit_length/2;
end    
for j=1:num_spatial
    %center position calculation
    x = 0; %will translate to fp axis system later
    y = -(num_spatial/2)*(d+d2) + (j-1)*(d+d2) + ((d+d2)/2);
    xy_s(1,j,1) = x;
    xy_s(2,j,1) = y;
    %corner position calculations
    for m=1:4
        xc = x + (slit_width/2) - slit_width*xyc(1,m); % alternates -d/2,+d/2,-d/2,+d/2
        yc = y + (d/2) - d*xyc(2,m); % alternates +d/2,+d/2,-d/2,-d/2
        xy_s(1,j,m+1) = xc;
        xy_s(2,j,m+1) = yc;
    end%for m
end%for j

%% definition of pixel positions, ideal initialization
xyz_p = zeros*[3,num_spectral,num_spatial,5]; %512 x,y,z coordinates for center, four corners
for j=1:num_spatial
    for k=1:num_spectral
        %center position calculation
        x = -(num_spectral/2)*d + (k-1)*d +(d/2);
        y = -(num_spatial/2)*(d+d2) + (j-1)*(d+d2) + ((d+d2)/2);
        z = 0;
        xyz_p(1,k,j,1) = x;
        xyz_p(2,k,j,1) = y;
        xyz_p(3,k,j,1) = z;
        %corner position calculations
        for m=1:4
            xc = x + (d/2)*xyc(1,m); % alternates -d/2,+d/2,-d/2,+d/2
            yc = y + (d/2)*xyc(2,m); % alternates +d/2,+d/2,-d/2,-d/2
            zc = 0.0;
            xyz_p(1,k,j,m+1) = xc;
            xyz_p(2,k,j,m+1) = yc;
            xyz_p(3,k,j,m+1) = zc;
        end%for m
    end %for k
end%for j


%% definition of "IFOV" pixel positions, ideal initialization
IFOVxyz_p = zeros*[3,num_spectral,num_spatial,5]; %512 x,y,z coordinates for center, four corners
for j=1:num_spatial
    for k=1:num_spectral
        %center position calculation
        x = -(num_spectral/2)*d + (k-1)*d +(d/2);
        y = -(num_spatial/2)*(d+d2) + (j-1)*(d+d2) + ((d+d2)/2);
        z = 0;
        IFOVxyz_p(1,k,j,1) = x;
        IFOVxyz_p(2,k,j,1) = y;
        IFOVxyz_p(3,k,j,1) = z;
        %corner position calculations
        for m=1:4
            xc = x + (slit_width/2)*xyc(1,m); % alternates -sw/2,+sw/2,-sw/2,+sw/2
            yc = y + (d/2)*xyc(2,m); % alternates +d/2,+d/2,-d/2,-d/2
            zc = 0.0;
            IFOVxyz_p(1,k,j,m+1) = xc;
            IFOVxyz_p(2,k,j,m+1) = yc;
            IFOVxyz_p(3,k,j,m+1) = zc;
        end%for m
    end %for k
end%for j


%% set up transformation for slit
rotation = [cos(tau_slit),sin(tau_slit);-sin(tau_slit),cos(tau_slit)];
translation = [x_slit, y_slit];
%% perform transformation for slit
for m=1:5
    xy = slit(:,m);
    xy = xy'*rotation+translation;
    slit(:,m) = xy;
end%for m

%% individual slit 'pixel's are not known really from slit transform, will project back positions from FP
%for j=1:num_spatial
%    for m=1:5
%        xy = xy_s(:,j,m);
%        xy = xy'*rotation+translation;
%        xy_s(:,j,m) = xy;
%    end%for m
%end%for j        

%testxs = reshape(xy_s(1,:,1),[8,1]);
%testys = reshape(xy_s(2,:,1),[8,1]);

%% set up transformation for fp
%rotation = [cos(tau_fp),sin(tau_fp),0;-sin(tau_fp),cos(tau_fp),0;0,0,1];
rotation_z = [cos(tau_fp),sin(tau_fp),0;-sin(tau_fp),cos(tau_fp),0;0,0,1];
rotation_y = [cos(tip_fp),0,sin(tip_fp);0,1,0;-sin(tip_fp),0,cos(tip_fp)];
rotation_x = [1,0,0;0,cos(tilt_fp),sin(tilt_fp);0,-sin(tilt_fp),cos(tilt_fp)];
rotation = rotation_x*rotation_y*rotation_z;
translation = [x_fp, y_fp, z_fp];
%% perform transformation for fp
for j=1:num_spatial
    for k=1:num_spectral
        for m=1:5
            xyz = xyz_p(:,k,j,m);
            %xyz = xyz'*rotation+translation;
            xyz = (xyz+translation')'*rotation; %expt with order of operations 10/17/22, this allows tip/tilt to be offcenter
            xyz_p(:,k,j,m) = xyz;

            xyz = IFOVxyz_p(:,k,j,m);
            xyz = (xyz+translation')'*rotation; %expt with order of operations 10/17/22, this allows tip/tilt to be offcenter
            IFOVxyz_p(:,k,j,m) = xyz;
        end%for m
    end %for k
end%for j        
%testxp = reshape(xyz_p(1,:,:,1),[num_spectral,num_spatial]);
%testyp = reshape(xyz_p(2,:,:,1),[num_spectral,num_spatial]);
%testzp = reshape(xyz_p(3,:,:,1),[num_spectral,num_spatial]);


%% project pixels back onto ideal focal (xy) plane with defocus
xy_p = zeros*[2,num_spectral,num_spatial,5]; %512 x,y coordinates for center, four corners
IFOVxy_p = zeros*[2,num_spectral,num_spatial,5]; %512 x,y coordinates for center, four corners
for j=1:num_spatial
    for k=1:num_spectral
        xy_p(1,k,j,1) = xyz_p(1,k,j,1);
        xy_p(2,k,j,1) = xyz_p(2,k,j,1);
        for m=1:4
            %defocus x and y anamorphically!
            x = xyz_p(1,k,j,m+1)+abs(xyz_p(3,k,j,m+1))*xyc(1,m)*tan(theta_spectral);
            y = xyz_p(2,k,j,m+1)+abs(xyz_p(3,k,j,m+1))*xyc(2,m)*tan(theta_spatial);
            xy_p(1,k,j,m+1) = x;
            xy_p(2,k,j,m+1) = y;
        end%for m

        IFOVxy_p(1,k,j,1) = IFOVxyz_p(1,k,j,1);
        IFOVxy_p(2,k,j,1) = IFOVxyz_p(2,k,j,1);
        for m=1:4
            %defocus x and y anamorphically!
            x = IFOVxyz_p(1,k,j,m+1)+abs(IFOVxyz_p(3,k,j,m+1))*xyc(1,m)*tan(theta_spectral);
            y = IFOVxyz_p(2,k,j,m+1)+abs(IFOVxyz_p(3,k,j,m+1))*xyc(2,m)*tan(theta_spatial);
            IFOVxy_p(1,k,j,m+1) = x;
            IFOVxy_p(2,k,j,m+1) = y;
        end%for m
    end %for k
end%for j

%scale image due to off-nominal focal lengths
for j=1:num_spatial
    for k=1:num_spectral
        for m=1:5
            xy_p(1,k,j,m) = TIRS.fl2*xy_p(1,k,j,m)/TIRS.fl2_nom;
            xy_p(2,k,j,m) = TIRS.fl*xy_p(2,k,j,m)/TIRS.fl_nom;

            IFOVxy_p(1,k,j,m) = TIRS.fl2*IFOVxy_p(1,k,j,m)/TIRS.fl2_nom;
            IFOVxy_p(2,k,j,m) = TIRS.fl*IFOVxy_p(2,k,j,m)/TIRS.fl_nom;
        end
    end
end

%% project each pixel through slit, determine relative offsets
%% want effective slit width and 'dispersion' distance, at each fp coordinate, also offset of projected image through the slit
%% perform projection of zero-order slit image onto FP
slit(1,:) = slit(1,:)+slit_image_position;
wavelengths = zeros([num_spectral, num_spatial,3]); %3 dimensions, for center, short edge, long edge
for j=1:num_spatial
    for k=1:num_spectral
        slit_here = slit(1,1)+xy_p(2,k,j,1)*tan(-tau_slit); %y position on fp determines clocking error in x, 
        dispersion_distance = xy_p(1,k,j,1)-slit_here;
        wavelengths(k,j,1) = dispersion_distance/dispersion;
        edge_left = (xy_p(1,k,j,2) + xy_p(1,k,j,4))/2; %left edges of each projected pixel averaged
        edge_right = (xy_p(1,k,j,3) + xy_p(1,k,j,5))/2; %right edges of each projected pixel averaged
        wavelengths(k,j,2) = (edge_left-slit_here)/dispersion;
        wavelengths(k,j,3) = (edge_right-slit_here)/dispersion;
    end %for k
end%for j        



%% routines for output
[outp_gm_fpath] = output_geometry_2(slit_width, d, num_spectral, ...
                   num_spatial, x_slit, y_slit, tau_slit, x_fp, y_fp, z_fp, ...
                   rad2deg(tau_fp), rad2deg(tip_fp), rad2deg(tilt_fp), xy_p, ...
                   slit, IFOVxy_p, TIRS_number);
[outp_wl_fpath] = output_wavelengths_2(slit_width, d, num_spectral, ...
                   num_spatial, x_slit, y_slit, tau_slit, x_fp, y_fp, z_fp, ...
                   rad2deg(tau_fp), rad2deg(tip_fp), rad2deg(tilt_fp), ...
                   d/dispersion, wavelengths, TIRS_number);
[gm] = load(outp_gm_fpath); %this just-created file is in 'the working' folder
[wl] = load(outp_wl_fpath); %this just-created file is in 'the working' folder
%% routines for plotting
if doplots
    plot_geometry(gm);
    plot_wavelengths(wl);
end

%% desire some statistics, like area of scene, mean and standard deviation
% until there is a projection onto earth, use 525 km
orbit_height = 525; % km
focal_length = TIRS.fl; % microns
area = zeros*[num_spectral,num_spatial]; %64x8
swap = [1 2 4 3]; % polyarea likes clockwise order
for k=1:num_spectral
    for j=1:num_spatial
        xpixel = reshape(xy_p(1,k,j,2:5),[1,4]); %use corners only
        ypixel = reshape(xy_p(2,k,j,2:5),[1,4]); %use corners only
        xpixel = xpixel(:,swap); %reorder to clockwise
        ypixel = ypixel(:,swap); %reorder to clockwise
        area(k,j) = polyarea(xpixel,ypixel);
    end
end

GSD = orbit_height*sqrt(area)/focal_length;
ratio = area/d^2;
max_GSD = max(max(GSD));
min_GSD = min(min(GSD));
avg_GSD = mean2(GSD);
std_GSD = std2(GSD);
inhomogeneity = std_GSD/avg_GSD;

end
