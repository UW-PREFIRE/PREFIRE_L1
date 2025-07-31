function [ground, geo_qf] = process_ground_footprints(TIRS, sat, E, geo_qf, ypr)
% ypr - yaw/pitch/roll angles, shaped (3,n)

%this function takes the precalculated azimuth and tilt angles and the
%spacecraft geodetic coordinates to determine the intersection points on
%the ellipsoid E, The az and tilt arguments specify the direction of the 
%view (line-of-sight) as the azimuth angle, measured clockwise from North, 
%and a tilt angle.
n_frames = length(sat.geod_lat);
ground.Ps = zeros([TIRS.spectral_channels,TIRS.spatial_scenes,5,n_frames]);
ground.Qs = zeros([TIRS.spectral_channels,TIRS.spatial_scenes,5,n_frames]);
ground.Pi = zeros([TIRS.spectral_channels,TIRS.spatial_scenes,5,n_frames]);
ground.Qi = zeros([TIRS.spectral_channels,TIRS.spatial_scenes,5,n_frames]);
% Note on array shapes: the lat/lon/alt are shaped (n,1); the az
% and tilt are (8,n); if we transpose the tilt and az, then it
% actually broadcasts as desired
tmp_nonneg_alt = max(sat.geod_alt, 0);
[ground.P,ground.Q,~] = lookAtSpheroid( ...
    sat.geod_lat, sat.geod_lon, tmp_nonneg_alt, ...
    sat.az', sat.tilt', E); %centers only

% For each frame, set boolean that determines which lat_ji_* to use for
%  trailing/leading edges.
%    Requires ypr(1,:) (i.e., yaw values) to be in the range of
%       (-90 to -270) or (90 to 270) degrees for flipped
%       (-90 to 90) degrees for normal
b_flipped = ((ypr(1,:) <= -90.) | (ypr(1,:) > 90.));

% the broadcasting will only work for 1 extra dimension (I think)
% so we need to loop over the outer two dimensions. Note that the
% tilt/az do not currently depend on the spectral channel, so the
% outermost loop is just replicating calculations at the moment.
tmp_nonneg_alt_b = max(sat.geod_alt_b, 0);
tmp_nonneg_alt_e = max(sat.geod_alt_e, 0);
for j = 1:TIRS.spectral_channels
    for i = 1:TIRS.spatial_scenes
        sqz_azs_T = squeeze(sat.azs(j,i,:,:))';
        sqz_tilts_T = squeeze(sat.tilts(j,i,:,:))';
        [lat_ji_b, lon_ji_b, ~] = lookAtSpheroid( ...
            sat.geod_lat_b, sat.geod_lon_b, tmp_nonneg_alt_b, ...
            sqz_azs_T, sqz_tilts_T, E);
        [lat_ji_c, lon_ji_c, ~] = lookAtSpheroid( ...
            sat.geod_lat, sat.geod_lon, tmp_nonneg_alt, ...
            sqz_azs_T, sqz_tilts_T, E);
        [lat_ji_e, lon_ji_e, ~] = lookAtSpheroid( ...
            sat.geod_lat_e, sat.geod_lon_e, tmp_nonneg_alt_e, ...
            sqz_azs_T, sqz_tilts_T, E);

        % For each frame, select which lat_ji_* to use for trailing/leading
        %  edges:
        lat_ji_tr = lat_ji_b;
        lat_ji_ld = lat_ji_e;
        lon_ji_tr = lon_ji_b;
        lon_ji_ld = lon_ji_e;
        lat_ji_tr(b_flipped,:) = lat_ji_e(b_flipped,:);
        lat_ji_ld(b_flipped,:) = lat_ji_b(b_flipped,:);
        lon_ji_tr(b_flipped,:) = lon_ji_e(b_flipped,:);
        lon_ji_ld(b_flipped,:) = lon_ji_b(b_flipped,:);

        lat_ji = lat_ji_tr;  % Init array, fill trailing edge vertices
        lat_ji(:,1) = lat_ji_c(:,1);   % Center
        lat_ji(:,2) = lat_ji_ld(:,2);  % Fill leading edge vertices
        lat_ji(:,4) = lat_ji_ld(:,4);  %
        lon_ji = lon_ji_tr;  % Init array, fill trailing edge vertices
        lon_ji(:,1) = lon_ji_c(:,1);   % Center
        lon_ji(:,2) = lon_ji_ld(:,2);  % Fill leading edge vertices
        lon_ji(:,4) = lon_ji_ld(:,4);  %

        ground.Ps(j,i,:,:) = lat_ji';
        ground.Qs(j,i,:,:) = lon_ji';

        % Save IFOV footprint corners at t = (TIRS image integration midpoint):
        ground.Pi(j,i,:,:) = lat_ji_c';
        ground.Qi(j,i,:,:) = lon_ji_c';
    end
end

end
