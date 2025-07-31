function [ground_corrected] = compute_view_angles( ...
    ground_corrected, sat, E_m, uc, tc)

% inputs:
% ground_corrected, a MATLAB structure various fields (the after
% terrain-correction positions), the following are used in this function:
%     lat: terrain-corrected scene latitude (nscene, nframe)
%     lon: terrain-corrected scene longitude (nscene, nframe)
%     alt: mean altitude of scene in meters. (nscene, nframe)
% sat, a MATLAB structure containing various fields related to the
% satellite position; the following are used in this function:
%     UTC_splitvals:  array of UTC time stamps, split into
%       [year, month, day, hour, minute, sec], shaped (nframe, 6)
%     ctimeN_minus_UTC: leap-second counts array, corresponds to (TAI-UTC) [s]
%     sat.geod_lat:  satellite geodetic latitude [deg_N], shape (nframe,1)
%     sat.geod_lon:  satellite geodetic longitude [deg_E], shape (nframe,1)
%     sat.geod_alt:  satellite geodetic altitude [km], shape (nframe,1)
% E_m, the referenceEllipsoid object, in units of meters.
% uc, the return from 'constants_unit_conversion'
% tc, the return from 'constants_time'
%
% outputs:
% various fields are added to the ground_corrected structure:
%     solar_zenith
%     solar_azimuth
%     solar_distance
%     sensor_zenith
%     sensor_azimuth
%%

lat = ground_corrected.lat;
lon = ground_corrected.lon;
alt = ground_corrected.alt;
[num_xtrack, num_atrack] = size(alt);

sat_lat = sat.geod_lat;
sat_lon = sat.geod_lon;
sat_alt = sat.geod_alt*uc.km_to_m;

solazi = zeros([num_xtrack, num_atrack]);
solzen = zeros([num_xtrack, num_atrack]);
soldst = zeros([num_xtrack, num_atrack]);

senazi = zeros([num_xtrack, num_atrack]);
senzen = zeros([num_xtrack, num_atrack]);

% eci2aer requires UTC input, along with the time deltas that are
% needed to convert to TT (Terrestrial Time). These are:
% deltaAT  = TAI - UTC
% deltaUT1 = UTC - UT1

% We don't know the deltaUT1 - ignoring it for now, since this is a
% < 1 second offset which means < 0.5 km geolocation error.
% deltaUT1 is the smooth variation of the time difference between
% UT1 and UTC, which is essentially the fractional second between
% TT and UT1.
deltaAT = sat.ctimeN_minus_UTC;
deltaUT1 = zeros([num_atrack, 1]);
% also assuming a zero polar motion, this could be retrieved from
% IERS data if we think the error due to ignoring it is "significant"
polarmotion = zeros([num_atrack,2]);

% for time processing, largely following this example:
%  https://www.mathworks.com/help/aerotbx/ug/
%      estimate-sun-analemma-using-planetary-ephemerides-and-eci-to-aer-transformation.html

% converting the UTC to approximate TAI (ignoring deltaUT1, only
% applying deltaAT, the leap seconds.)
UTC_dts = datetime(sat.UTC_splitvals);
approx_TAI_dts = UTC_dts + seconds(deltaAT);

% converting approx TAI to TT
approx_TT_dts = approx_TAI_dts + seconds(tc.TT_minus_TAI);

% convert the TT (as datetime objects) to the julian date TDB
% needed by planetEphemeris, which returns ECI position of the sun
% note we convert match into a 'date vector' here, as that is what
% the builtin MATLAB function wants.
approx_TT_datevec = [year(approx_TT_dts), ...
                    month(approx_TT_dts), ...
                    day(approx_TT_dts), ...
                    hour(approx_TT_dts), ...
                    minute(approx_TT_dts), ...
                    second(approx_TT_dts)];
TDB_jdates = tdbjuliandate(approx_TT_datevec);

% choice of Ephemeris model (405) is the default, and defines
% positions with respect to ICRS 1. See:
% https://www.mathworks.com/help/aerotbx/ug/planetephemeris.html
sun_eci = planetEphemeris(TDB_jdates, 'Earth', 'Sun', 405, 'km');
sun_eci = sun_eci * uc.km_to_m;

% This could be done with a single flattened array, but the
% implementation is simpler to do as a loop over the xtrack positions.
% This loop reuses the solar position, which is constant for all
% xtrack positions at one atrack time.
for x=1:num_xtrack

    % use eci2aer for solar position (since we get the Sun position
    % as an ECI coordinate from the planetEphemeris);
    % geodetic2aer for the sensor position (since the satellite
    % position is in geodetic lat/lon/alt.)
    %
    % these come from different MATLAB toolboxes so the
    % input/output array expectations are different: eci2aer wants
    % a (n,3) computed array of lla, while geodetic2aer wants
    % separate arrays.
    lla = [lat(x,:); lon(x,:); alt(x,:)]';
    aer = eci2aer( ...
        sun_eci, sat.UTC_splitvals, lla, ...
        'IAU-2000/2006', deltaAT, deltaUT1, polarmotion, ...
        'flattening', E_m.Flattening, ...
        're', E_m.SemimajorAxis);
    solazi(x,:) = aer(:,1)';
    solzen(x,:) = 90 - aer(:,2)';
    soldst(x,:) = aer(:,3)';  % [m]

    % not saving the slantrange.
    [azim, elev, ~] = geodetic2aer( ...
        sat.geod_lat, sat.geod_lon, sat_alt, ...
        lla(:,1), lla(:,2), lla(:,3), ...
        E_m);
    senazi(x,:) = azim;
    senzen(x,:) = 90 - elev;

end

ground_corrected.solar_azimuth  = solazi;
ground_corrected.solar_zenith   = solzen;
ground_corrected.solar_distance = soldst;  % [m]

ground_corrected.sensor_azimuth = senazi;
ground_corrected.sensor_zenith  = senzen;

end
