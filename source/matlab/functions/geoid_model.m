%
% Geoid implementation
%
% This is implemented as an object with a similar layout to the
% DEM, in order to make the usage of the geoid and DEM behave
% similarly. The geoid does not really need this approach, since
% the whole thing lives in one file, but the class is then very simple.
%
%
classdef geoid_model < handle

properties(SetAccess = private)
    geoid_file
    z
    info
    lon0
    lat0
    dlon
    dlat
    nrow
    ncol
    units
end

methods

    function obj = geoid_model(geoid_geotiff)
    %geoid_model Object class to facilitate access to Geoid model.
    %
    % geoid_obj = geoid_model(geoid_geotiff)
    % creates the object given an input geotiff file, containing
    % the geoid model data.

        obj.load_geoid(geoid_geotiff);
    end

    function geoid_z = lookup(obj, lon, lat)
    % geoid_model.lookup find the z values according to input lon/lat values
    %
    % geoid_z = geoid_model.lookup(lon, lat)
    % returns z values, in meters, from the geoid surface sampled
    % at the input lon, lat locations. lon & lat can be arrays of
    % any (matching) sizes, in which case the returned geoid_z
    % will have the same size. The lookup is done by fast integer
    % lookup: divide by increment, then floor.

        % note that the geoid data array is stored with latitude the
        % first axis (rows) and longitude the second (columns)
        r = floor((lat - obj.lat0) / obj.dlat) + 1;
        c = floor((lon - obj.lon0) / obj.dlon) + 1;

        rmsk = (r > obj.nrow) | (r < 1);
        cmsk = (c > obj.ncol) | (c < 1);

        if any(rmsk)
            r(rmsk) = 1;
        end
        if any(cmsk)
            c(cmsk) = 1;
        end

        if any(rmsk) || any(cmsk)
            fprintf('Warning, out of range lookup positions set to nan in output\n');
        end

        ind = sub2ind(size(obj.z), r, c);
        geoid_z = obj.z(ind);

        if any(rmsk)
            geoid_z(rmsk) = nan;
        end
        if any(cmsk)
            geoid_z(cmsk) = nan;
        end

    end

end

methods(Access = private)

    function load_geoid(obj, geoid_geotiff)
    % helper function to read the geotiff file, and store the data
    % in various object properties.
        obj.geoid_file = geoid_geotiff;
        % note the file is assumed to contain the geoid in units of
        % meters. True for "us_nga*.tiff" files from Agisoft.
        obj.z = double(imread(geoid_geotiff));
        obj.info = geotiffinfo(geoid_geotiff);
        % copy out some fields from geotiff info into direct object
        % attributes for ease of use.
        % I don't think this is the "correct, general approach", but this
        % was where I could find the numbers that made sense for the
        % lat/lon grids.
        obj.lon0 = obj.info.GeoTIFFTags.ModelTiepointTag(4);
        obj.lat0 = obj.info.GeoTIFFTags.ModelTiepointTag(5);
        obj.dlon =  obj.info.PixelScale(1);
        obj.dlat = -obj.info.PixelScale(2);
        obj.nrow = size(obj.z, 1);
        obj.ncol = size(obj.z, 2);
        obj.units = 'meters';
    end

end

end
