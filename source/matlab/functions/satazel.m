function [sat] = satazel(sat, ypr, TIRS, gm)
%
% calculate the azimuth and tilt angles for all spatial/spectral scenes.
%
% INPUT:
% sat is a MATLAB structure, that must contain one field:
% beta, an angle in degrees, shaped (1, n) for n frames.
% beta is not clearly defined anywhere, so the exact definition is TBD.
%
% ypr - yaw/pitch/roll angles, shaped (3,n). These are applied to
% the attitude data in sat, and the instrument geometry.
%
% TIRS - the MATLAB structure containing various instrument parameters
% gm - the MATLAB structure loaded from the static geometry file
%
% OUTPUT:
% the following fields are added to the sat structure:
% tilts - tilt angles in degrees, for the footprint centers and
%     corner vertices. shape (64,8,5,n)
% azs - azimuth angles in degrees, for the footprint centers and
%     corner vertices. shape (64,8,5,n)
% tilt - tilt angle for the footprint center, first channel, shape (8,n)
% az - az angle for the footprint center, first channel, shape (8,n)
%
% Note: "tilt" and "azimuth" are direct inputs to the MATLAB
% lookAtSpheroid function.
% These describe the line of sight vector from the observer at
% lat/lon/alt; the "tilt" angle is what I think would be normally
% called the nadir angle (angle between LOS and nadir); the azimuth
% is the normally defined azimuth (clockwise from 0 deg North, as
% viewed from outside the Earth looking down), of the projection of
% the LOS down to the Earth surface tangent at the observer
% lat/lon. See:
% https://www.mathworks.com/help/map/ref/lookatspheroid.html

%(TIM) * 'yaw' should be a rotation about the CubeSat y-axis (which is oriented
%         normal to and toward the exterior of the "nadir plate")
%      * 'pitch' should be a rotation about the CubeSat z-axis. A communication-
%         session reorientation is a rotation about the CubeSat z-axis.
%      * 'roll' should be a rotation about the CubeSat negative x-axis (i.e.,
%         the science-mode along-orbit direction of movement)

%pointing for center and 4 corners/center from geometry file
%recall gm.IFOVxy_p are coordinates of each pixel in microns, 2 x 64 x 8 x 5 = xy * spectral * spatial * center & 4 corners
%need an 'ideal' metric for translation of x (spectral) pixel positions
xd = (0:TIRS.spectral_channels-1)*TIRS.d+TIRS.slit_image_position; 
num_frames = size(sat.beta,2);
% get an attitude adjustment and calculate look angles / would be great to vectorize this 
for i=1:TIRS.spatial_scenes
    for j=1:TIRS.spectral_channels
        for k=1:5 %center and four corners
            alphas(j,i,k,:) = deg2rad(ypr(3,:))+(gm.IFOVxy_p(2,j,i,k))/TIRS.fl; %add roll
            for m = 1:num_frames
                rel_x(j,i,k,m) =  deg2rad(ypr(2,m))+TIRS.smile(i)+((gm.IFOVxy_p(1,j,i,k) - xd(j))/TIRS.fl2); %should be place to add pitch, anamorphic optics require two different focal lengths
            end
            alpha_primes(j,i,k,:) = sign(alphas(j,i,k,:)).*rad2deg(sqrt(alphas(j,i,k,:).*alphas(j,i,k,:) + rel_x(j,i,k,:).*rel_x(j,i,k,:)));
            beta_primes(j,i,k,:) = -rad2deg(atan(rel_x(j,i,k,:)./alphas(j,i,k,:))); 
            leftrights(j,i,k,:) = -0.5*(1+sign(alphas(j,i,k,:)))*180;%[0 0 0 0 -180 -180 -180 -180]; %swap angle of azimuth for each side of swath, could be wrong if roll is greater than alpha
            sat.azs(j,i,k,:) = sat.beta(1,:) + reshape(leftrights(j,i,k,:), [1,num_frames]) + ypr(1,:) + reshape(beta_primes(j,i,k,:), [1,num_frames]); % beta may already contain yaw, if so, it is double counted
            sat.tilts(j,i,k,:) = abs(reshape(alpha_primes(j,i,k,:),[1,num_frames])); %pitch inside of absolute value *should* work
        end
    end
end

%centers only
sat.az = reshape(sat.azs(1,:,1,:),[TIRS.spatial_scenes,num_frames]);
sat.tilt = reshape(sat.tilts(1,:,1,:),[TIRS.spatial_scenes,num_frames]);

end
