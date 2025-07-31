%
% Elevation model implementation
%
% general idea here is to implement an object where a caller can
% just request dem data at a certain location, either by nearest
% neighbor lookup (for a single value) or to return an array of DEM
% elements (for a lat/lon polygon). The object class internally
% manages access to the DEM tiles, so that the caller does not need
% to know any of those details.
%
% We use an optional sized buffer of tiles, where requests to lat/lon
% positions that do not fall in any currently loaded tile will
% cause a new tile to be loaded.
% If the buffer is already full, then the "most stale" tile is
% deleted, to conserve memory. "Most stale" means that the last
% access of that tile was most separated from the current access -
% in other words, that tile has not been used in the longest time.
% This also means that the method will be most efficient with a
% buffer size of at least 4 (so that lookups at edges or corners
% will be efficiently handled), and that the loop over lat/lon
% positions that the caller is doing should be ordered in some way
% that nearby lookups are done in successive iterations. For
% example, if the lookups are following a satellite ground track,
% then just order the calls according to the along track distance
% (or time).
% For wide - swath data, the number of buffered tiles should be
% increased, so that the entire cross-track extent at one time can
% be covered by the loaded DEM tiles.
%
%
% this manages both elevation tiles and land-water mask. (TBD: height error?)
% designed to work with copernicus DEM @ 90 m.
%
% WBM mask from copernicus
% 0 : no water (land)
% 1 : ocean
% 2 : lake
% 3 : river
%

classdef elevation_model < handle

% note: move properties to hidden if they are not generally useful
% for interactive work or testing.
% Note that hidden variables can still be accessed at the console,
% but cannot be 'discovered' by tab completion.
properties (SetAccess = private)
    root_dir
    buffer_size
    max_buffer_size
    tile_ids
    tile_elev
    tile_info
    tile_mask
    tile_is_ocean
    tile_ll
    % currently, the tile coverage (the corners) is not
    % used. Leaving it here, but commented, in case this is
    % needed again.
    %tile_coverage_ul
    %tile_coverage_lr
    tile_polyshapes
    tile_efilepaths
    tile_mfilepaths
    last_access
    write_access_log
    access_log_file
    verbose
end

properties (SetAccess = private, Hidden)
    debug
end

methods

    function obj = elevation_model( ...
        root_dir, max_buffer_size, ...
        verbose, debug, access_log_file)
    %elevation_model Object class to facilitate access to a tiled Digital Elevation Model (DEM)
    %
    % dem_obj = elevation_model(root_dir, max_buffer_size)
    % creates the object with a given buffer size (in number of
    % tiles), using the root_dir as the location of the DEM tile
    % data. Thie buffer_size limits the number of tiles that will be loaded
    % into memory at once.
    %
    %
    % dem_obj = elevation_model(root_dir, max_buffer_size, verbose)
    % can be used to turn on console printed information (via a
    % logical value true) when a tile is loaded, and what the
    % average access count is. This is helpful to assess performance
    % with target testing or iteractive work.
    %
    % dem_obj = elevation_model(root_dir, max_buffer_size, verbose, debug)
    % can be used to turn on a debug mode (via a logical value
    % true) that prints information to console, primarily to track
    % the order of execution of methods.
    %
    % dem_obj = elevation_model(root_dir, max_buffer_size, verbose, debug, log_file)
    % the final optional input, log_file, is a string file+path for
    % a log file, that contains an updated line of the buffer
    % contents each time a new tile is loaded. Each line contains
    % the buffer_id and the last_access counter. This log file can
    % be used to analyze the tile usage pattern.
        if nargin >= 3
            obj.verbose = verbose;
        else
            obj.verbose = false;
        end
        if nargin >= 4
            obj.debug = debug;
        else
            obj.debug = false;
        end
        if nargin == 5
            obj.write_access_log = true;
            obj.access_log_file = access_log_file;
        else
            obj.write_access_log = false;
        end

        if obj.debug
            disp('->elevation_model');
        end

        % for copernicus DEM
        obj.root_dir = root_dir;

        obj.buffer_size = 0;
        obj.max_buffer_size = max_buffer_size;
        % in MATLAB, row vectors seem to behave most like 1D vectors.
        obj.tile_ids = cell(1,max_buffer_size);
        obj.tile_elev = cell(1,max_buffer_size);
        obj.tile_info = cell(1,max_buffer_size);
        obj.tile_mask = cell(1,max_buffer_size);
        obj.tile_is_ocean = cell(1,max_buffer_size);
        obj.tile_ll = cell(1,max_buffer_size);
        % see note in property definitions about tile_coverage.
        %obj.tile_coverage_ul = zeros(2,max_buffer_size);
        %obj.tile_coverage_lr = zeros(2,max_buffer_size);
        obj.tile_polyshapes = cell(1,max_buffer_size);
        obj.tile_efilepaths = cell(1,max_buffer_size);
        obj.tile_mfilepaths = cell(1,max_buffer_size);
        obj.last_access = zeros(1,max_buffer_size);

        for n=1:max_buffer_size
            obj.tile_ids{n} = 'empty';
        end

    end

    function [h, m] = get_nearest_neighbor_h(obj, lon, lat)
    %get_nearest_neighbor_h simplest DEM lookup, which finds the h and m values for the nearest lon/lat value.
    %
    % [h, m] = dem.get_nearest_neighbor_h(lon, lat)
    % for a given lat, lon position (in degrees), returns the
    % surface altitude in meters (h) and land/water mask value (m)
    % for the nearest point in the DEM.
        if obj.debug
            disp('->get_nearest_neighbor_h');
        end

        slot = obj.find_buffer_by_ll(lon, lat);
        if obj.tile_is_ocean{slot}
            h = 0.0;
            m = 0;
        else
            h = elevation_model.nearest_neighbor_h( ...
                obj.tile_elev{slot}, ...
                obj.tile_info{slot}, ...
                lon, lat);
            m = elevation_model.nearest_neighbor_h( ...
                obj.tile_mask{slot}, ...
                obj.tile_info{slot}, ...
                lon, lat);
        end
    end

    function [h, m, wt, ct] = get_pixel_group_h(obj, P)
    %get_pixel_group_h get alt and mask values for all DEM pixels inside polyshape P
    %
    % returns two 1D arrays of the same shape, containing the
    % altitudes (h) and the land / water mask (m) for all DEM
    % pixels inside P.
    % third 1D array (wt) contains the weights; computing the mean
    % value of h or m can be done as sum(h * w) instead of mean(h),
    % in order to properly account for DEM resolution changes at
    % tile edges.
    %
    % the first three arrays are all the same shape, containing (k)
    % values from the high resolution DEM data.
    %
    % the fourth (and optional) array ct contains the number of
    % overlapping DEM points from each overlapping tiles. In this
    % case, for efficiency, only one array of the point counts from
    % each of the N tiles is returned, meaning a 1D array
    % containing (t) values where t is the number of tiles,
    % typically of order 1 - 10. This output is mainly for
    % debugging, most calls won't need these values.

        if obj.debug
            disp('->get_pixel_group_h');
        end

        [Plist, Pselection, slots] = obj.resolve_tile_matchups(P);

        n_slots = length(slots);
        hval_s = cell(n_slots,1);
        mval_s = cell(n_slots,1);
        wval_s = cell(n_slots,1);
        ct = zeros(n_slots,1);

        for n = 1:n_slots
            p = Pselection(n);
            [hval_s{n}, mval_s{n}, Pfrac] = ...
                obj.get_pixel_group_h_direct(Plist{p}, slots(n));
            ct(n) = numel(hval_s{n});
            wval_s{n} = ones(size(hval_s{n})) * Pfrac / ct(n);
        end

        if n_slots > 1
            h = cell2mat(hval_s);
            m = cell2mat(mval_s);
            wt = cell2mat(wval_s);
        else
            h = hval_s{1};
            m = mval_s{1};
            wt = wval_s{1};
        end

    end

    function [lon, lat] = get_pixel_group_ll(obj, P)
    %get_pixel_group_h get lon and lat for all DEM pixels inside polyshape P
    %
    % returns two 1D arrays of the same shape, containing the
    % longitude and latitude for all DEM pixels inside P.
    %
    % This is implemented as a separate method, because in
    % operations this is not currently needed. At the moment the
    % DEM pixel lat/lon would only be used for detailed
    % visualization of the interior altitude or mask of the polyshape.
    %
    % the two returned arrays will be of the same size as the
    % arrays returned from get_pixel_group_h.
        if obj.debug
            disp('->get_pixel_group_ll');
        end

        [Plist, Pselection, slots] = obj.resolve_tile_matchups(P);

        n_slots = length(slots);
        lon_s = cell(n_slots,1);
        lat_s = cell(n_slots,1);

        for n = 1:n_slots
            p = Pselection(n);
            % we can discard the Pfrac here, not needed for weighting
            [lon_s{n}, lat_s{n}, ~] = ...
                obj.get_pixel_group_ll_direct(Plist{p}, slots(n));
        end

        if n_slots > 1
            lon = cell2mat(lon_s);
            lat = cell2mat(lat_s);
        else
            lon = lon_s{1};
            lat = lat_s{1};
        end

    end

end

methods (Access = private)

    function [Plist, Pselection, slots] = resolve_tile_matchups(obj, P)
    %resolve_tile_matchups locates the DEM tiles that are needed to overlap with the input polyshape P
    %
    % This helper function contains various function calls to
    % manage corner cases dealiing with tile boundaries and the
    % antimeridian. To handle the antimeridian, the polyshape might
    % be "twinned" (to make an East and West polygon)
    % the returned arrays are the list of Polygons (either 1 or 2,
    % if the polygon was twinned at the antimeridian), an array of
    % selections from the polygon list, and an array of slot
    % numbers. These arrays are organized such that the pixel
    % lookup based on slots(n) should use the polygon
    % Plist(Pselection(n)).
    % Plist is a cell array of polyshapes; Pselection and slots are
    % both integer arrays.
    %
    % This function will also cause new tiles to be loaded as
    % needed, in order to fully cover the polyshape P; note that
    % the method will fail if P is so large that is requires a
    % number of tiles larger than the buffer size.

        antimerid = elevation_model.contains_antimeridian(P);

        if antimerid
            [Eslots, Wslots] = obj.find_slots_by_polyshape_antimerid(P);
            [EP, WP] = elevation_model.poly_antimeridian_unfold(P);
            slots = [Eslots; Wslots];
            Pselection = ones(length(Eslots)+length(Wslots),1);
            Pselection(length(Eslots)+1:end) = 2;
            Plist = {EP, WP};
        else
            slots = obj.find_slots_by_polyshape_direct(P);
            Pselection = ones(length(slots),1);
            Plist = {P};
        end

    end

    function slots = find_slots_by_polyshape(obj, P)
    %find_slots_by_polyshape find tile slots for an input polyshape
    %
    % this is mainly a convenience function for testing. in
    % operational processing, this step will be performed by the
    % resolve_tile_matchups() method, and it will call the direct
    % or antimerid versions as appropriate.
    %
    % Note, this function will return a list of slots needed to cover
    % the polyshape, and will trigger tile loads as needed.
        if elevation_model.contains_antimeridian(P)
            [Eslots, Wslots] = obj.find_slots_by_polyshape_antimerid(P);
            slots = [Eslots; Wslots];
        else
            slots = obj.find_slots_by_polyshape_direct(P);
        end
    end

    function slots = find_slots_by_polyshape_direct(obj, P)
    %find_slots_by_polyshape_direct find tile slots for an input polyshape
    %
    % Not intended to be used externally. This function assumes the
    % polyshape does not cross the antimeridian - it will give
    % wrong results in that case.
    % returns an array of slot numbers, and new tiles are loaded if needed.
        if obj.debug
            disp('->find_slots_by_polyshape_direct_');
        end

        lon_min = min(P.Vertices(:,1));
        lon_max = max(P.Vertices(:,1));
        lat_min = min(P.Vertices(:,2));
        lat_max = max(P.Vertices(:,2));
        [tlon1, tlat1] = elevation_model.ll_to_copernicus_tilell( ...
            lon_min, lat_min, 0.0);
        [tlon2, tlat2] = elevation_model.ll_to_copernicus_tilell( ...
            lon_max, lat_max, 0.0);
        s = 1;
        slots = zeros( (tlon2-tlon1+1)*(tlat2-tlat1+1), 1 );

        if obj.debug
            disp('direct tile lon1 lon2 lat1 lat2:')
            disp([tlon1 tlon2 tlat1 tlat2]);
        end

        for tlon = tlon1:tlon2
            for tlat = tlat1:tlat2
                tile_id = elevation_model.tilell_to_id(tlon, tlat);
                slots(s) = obj.find_buffer_by_tile_id(tile_id);
                s = s + 1;
            end
        end
    end

    function [east_slots, west_slots] = find_slots_by_polyshape_antimerid(obj, P)
    %find_slots_by_polyshape_antimerid find tile slots for an input polyshape
    %
    % Not intended to be used externally. This function assumes the
    % polyshape crosses the antimeridian - it will give wrong
    % results if the polyshape does not cross the antimeridian.
    % returns two arrays of slot numbers, for the requires East
    % tiles (near E179) and West tiles (W180, etc).
    % New tiles are loaded if needed.
        if obj.debug
            disp('->find_slots_by_polyshape_antimerid_');
        end

        lon = P.Vertices(:,1);
        east_lon = lon(lon >= 0);
        west_lon = lon(lon < 0);
        east_lon_min = min(east_lon);
        west_lon_max = max(west_lon);
        lat_min = min(P.Vertices(:,2));
        lat_max = max(P.Vertices(:,2));
        [tlon1, tlat1] = elevation_model.ll_to_copernicus_tilell( ...
            east_lon_min, lat_min, 0.0);
        % dummy point to basically find the E179 tile
        [tlon2, ~] = elevation_model.ll_to_copernicus_tilell( ...
            179.9, lat_min, 0.0);
        % dummy point to basically find the W180 tile
        [tlon3, ~] = elevation_model.ll_to_copernicus_tilell( ...
            -179.9, lat_max, 0.0);
        [tlon4, tlat2] = elevation_model.ll_to_copernicus_tilell( ...
            west_lon_max, lat_max, 0.0);

        east_slots = zeros( (tlon2-tlon1+1)*(tlat2-tlat1+1), 1 );
        west_slots = zeros( (tlon4-tlon3+1)*(tlat2-tlat1+1), 1 );
        if obj.debug
            disp('antimerid tile lon1 lon2 lon3 lon4 lat1 lat2:')
            disp([tlon1 tlon2 tlon3 tlon4 tlat1 tlat2]);
        end

        s = 1;
        for tlon = tlon1:tlon2
            for tlat = tlat1:tlat2
                tile_id = elevation_model.tilell_to_id(tlon, tlat);
                east_slots(s) = obj.find_buffer_by_tile_id(tile_id);
                s = s + 1;
            end
        end

        s = 1;
        for tlon = tlon3:tlon4
            for tlat = tlat1:tlat2
                tile_id = elevation_model.tilell_to_id(tlon, tlat);
                west_slots(s) = obj.find_buffer_by_tile_id(tile_id);
                s = s + 1;
            end
        end
    end


    function [h, m, Pfrac] = get_pixel_group_h_direct(obj, P, slot)
    %get_pixel_group_h_direct new core function for alt and mask lookup
    %
    % Not intended to be called externally. In this method, we assume
    % that the input polyshape does not need to be split, but can
    % always by intersected directly to the given tile (specified
    % via slot number.)
    % returns an array of altitudes and mask values (h and m),
    % and a scalar Pfrac, which is the fraction of P that resides
    % within the tile at the given slot number.
        Pn = intersect(P, obj.tile_polyshapes{slot});
        Pfrac = Pn.area/P.area;
        % here we intersect with the original P. Does it need to be Pn?
        % I think that should be the same result, but perhaps P is simpler,
        % since it would have 4 sides, but Pn could have 4-6 sides?
        if obj.tile_is_ocean{slot}
            h = elevation_model.extract_pixel_group( ...
                obj.tile_elev{slot}, ...
                obj.tile_ll{slot}, ...
                obj.tile_info{slot}, ...
                P);
            m = zeros(size(h), 'logical');
        else
            [h, m] = elevation_model.extract_pixel_group( ...
                obj.tile_elev{slot}, ...
                obj.tile_ll{slot}, ...
                obj.tile_info{slot}, ...
                P, ...
                obj.tile_mask{slot});
        end
    end

    function [lon, lat, Pfrac] = get_pixel_group_ll_direct(obj, P, slot)
    %get_pixel_group_ll_direct new core function for lon/lat lookup
    %
    % Not intended to be called directly. In this method, we assume
    % that the input polyshape does not need to be split, but can
    % always by intersected directly to the given tile (specified
    % via slot number.)
    % returns an array of longitudes and latitudes.
    % and a scalar Pfrac, which is the fraction of P that resides
    % within the tile at the given slot number.
        Pn = intersect(P, obj.tile_polyshapes{slot});
        Pfrac = Pn.area/P.area;
        [lon, lat] = elevation_model.extract_pixel_group( ...
            obj.tile_ll{slot}(:,:,1), ...
            obj.tile_ll{slot}, ...
            obj.tile_info{slot}, ...
            P, ...
            obj.tile_ll{slot}(:,:,2));
    end

    function slot = find_buffer_by_ll(obj, lon, lat)
    %find_buffer_by_ll return the buffer slot for the tile that contains the input lon,lat
    %
    % the tiles are basically tracked by converting the input
    % lon/lat (an arbitrary value) to the tile lon/lat (quantized
    % by 1x1 degree tiles), and then into a tile_id, a string of
    % the form "N12_E123". This is then passed to
    % find_buffer_by_tile_id.
    % returns the scalar integer slot number for the tile.
        if obj.debug
            disp(['->find_buffer_by_ll, lon: ', num2str(lon), ...
                  ' lat: ', num2str(lat)]);
        end

        % In order to know the tile lon/lat, (and in turn, the id)
        % we need to assume a grid. To do this, we call a function
        % that has a hardcoded knowledge of the copernicus grid
        [tlon, tlat] = elevation_model.ll_to_copernicus_tilell(lon,lat);
        tile_id = elevation_model.tilell_to_id(tlon, tlat);

        slot = obj.find_buffer_by_tile_id(tile_id);

    end

    function slot = find_buffer_by_tile_id(obj, tile_id)
    %find_buffer_by_tile_id return the buffer slot for the input tile_id
    %
    % if the requested tile_id is in the buffer already, then its
    % slot number is returned. If the requested tile_id is not in
    % the buffer, it is loaded by calling load_tile.
    % and the end, the access counter is incremented up by one for
    % all other buffered tiles.
        if obj.debug
            disp(['->find_buffer_by_tile_id: ', tile_id]);
        end

        slot = -1;
        for n=1:obj.buffer_size
            if ( obj.tile_ids{n} == tile_id )
                slot = n;
                break
            end
        end
        if slot == -1
            % since a tile containing that lon/lat is not loaded,
            % we need to load a new buffer.
            slot = obj.load_tile(tile_id);
        end
        % increment access counter by 1, for all tiles other than
        % the one used for this request.
        % this method could be sub-optimal, if the find_buffer
        % method is called often in a way that is not related to
        % pixel lookups. this is not currently the case.
        for s=1:obj.buffer_size
            if s == slot
                obj.last_access(s) = 0;
            else
                obj.last_access(s) = obj.last_access(s) + 1;
            end
        end

    end

    function slot = find_new_bufferslot(obj)
    %find_new_bufferslot find the buffer slot for a new tile.
    %
    % if the buffer is not yet full, this will return a new slot,
    % which implies the buffer will be extended; otherwise it
    % returns the slot that has the oldest last access.
        if obj.debug
            disp('-> find_new_bufferslot');
        end
        % if the buffer is not full, pick a new slot at the end of
        % the buffer cell arrays.
        % if the buffer is full, the new slot is selected as the
        % one with the oldest last access.
        if obj.buffer_size < obj.max_buffer_size
            obj.buffer_size = obj.buffer_size + 1;
            slot = obj.buffer_size;
        else
            [~,slot] = max(obj.last_access);
        end
    end

    function slot = load_tile(obj, tile_id)
    %load_tile load a new tile, given an input tile_id, and return the slot used for the new tile.
    %
    % find_new_bufferslot is used to get a new slot position for
    % this tile.
    % returns the integer slot number that was used for the loaded
    % tile.
    % several helper functions are called to handle the tile
    % data. If no tile file exists for this tile_id, it is assumed
    % to be an ocean tile, and a 'dummy' zero altitude tile is
    % loaded in place.
        if obj.debug
            disp(['->load_tile, tile_id: ' tile_id])
        end

        [efilepath, mfilepath] = ...
            elevation_model.construct_filepaths(obj.root_dir, tile_id);

        if exist(efilepath, 'file') == 2
            if obj.debug
                disp(['  reading tile: ' efilepath]);
            end
            [elev,info] = readgeoraster(efilepath);
            [wbm_mask,~] = readgeoraster(mfilepath);
            is_ocean = false;
            % land is integer 0, other values (1, 2, 3) are various
            % water bodies, (ocean, lake, river).
            mask = wbm_mask == 0;
        else
            if obj.debug
                disp(['  dummy tile for id: ' tile_id]);
            end
            % if the file does not exist, then create a fake 'ocean tile'
            % that contains zero terrain elevation, but with the same
            % array dimensions. The bounding corners are created for
            % what would be the tile in this position.
            info = elevation_model.make_dummy_tile_info(tile_id);
            elev = zeros(info.RasterSize,'single');
            mask = zeros(info.RasterSize,'logical');
            is_ocean = true;
            efilepath = 'ocean';
            mfilepath = 'ocean';
        end

        % for initial testing, always print when a tile was loaded.
        % this should be removed later when the code seems more stable.
        avg_last_access = mean(obj.last_access);
        if obj.verbose
            if is_ocean
                disp(['-> load_tile, ' ...
                      sprintf('avg last access %7.0f, ', avg_last_access) ...
                      'tile_id: ' tile_id ' (ocean tile)'])
            else
                disp(['-> load_tile, ' ...
                      sprintf('avg last access %7.0f, ', avg_last_access) ...
                      'tile_id: ' tile_id])
            end
        end

        slot = obj.find_new_bufferslot();
        obj.tile_ids{slot} = tile_id;
        obj.tile_elev{slot} = elev;
        obj.tile_mask{slot} = mask;
        obj.tile_info{slot} = info;
        obj.tile_ll{slot} = elevation_model.make_ll_grid(info);
        obj.tile_efilepaths{slot} = efilepath;
        obj.tile_mfilepaths{slot} = mfilepath;
        obj.last_access(slot) = 0;
        obj.tile_is_ocean{slot} = is_ocean;
        % see note in property definitions about tile_coverage.
        %[tile_ul, tile_lr] = elevation_model.compute_tile_coverage(info);
        %obj.tile_coverage_ul(:,slot) = tile_ul;
        %obj.tile_coverage_lr(:,slot) = tile_lr;
        obj.tile_polyshapes{slot} = elevation_model.compute_tile_polyshape(info);

        if obj.write_access_log
            fid = fopen(obj.access_log_file, 'a');
            for n=1:obj.max_buffer_size-1
                fprintf(fid, '%8s, %8d,', obj.tile_ids{n}, obj.last_access(n));
            end
            n = obj.max_buffer_size;
            fprintf(fid, '%8s, %8d\n', obj.tile_ids{n}, obj.last_access(n));
            fclose(fid);
        end
    end

end

% static methods deal with the actual extract of DEM points from
% the tiles, and other misc calcs, that do not need internal object data.
methods(Static)

    function [tfilepath, mfilepath] = construct_filepaths( ...
        root_dir, tile_id)
    % construct_filepaths construct filepath for the tile that resides in root_dir, with input tile_id
    %
    % returns both the terrain file (tfilepath) and the mask file
    % (mfilepath). Note that the names are created based on the
    % tile_id, which is set by some lon/lat location. If the
    % location is in the ocean, the tile may not exist.

        % examples
        % Copernicus_DSM_COG_30_N30_00_E030_00_DEM
        % Copernicus_DSM_COG_30_N30_00_E030_00_DEM/ \
        %     Copernicus_DSM_COG_30_N30_00_E030_00_DEM.tif
        % Copernicus_DSM_COG_30_N30_00_E030_00_DEM/AUXFILES \
        %     Copernicus_DSM_COG_30_N30_00_E030_00_WBM.tif
        %
        % in these cases, tile_id == N30_E030
        latpart = tile_id(1:3);
        lonpart = tile_id(5:8);
        nameroot = ['Copernicus_DSM_COG_30_' ...
                    latpart '_00_' lonpart '_00'];

        dem_subdir = [nameroot '_DEM'];
        tfilename = [nameroot '_DEM.tif'];
        tfilepath = [root_dir filesep dem_subdir filesep tfilename];

        wbm_subdir = [nameroot '_DEM' filesep 'AUXFILES'];
        mfilename = [nameroot '_WBM.tif'];
        mfilepath = [root_dir filesep wbm_subdir filesep mfilename];

    end

    function info = make_dummy_tile_info(tile_id)
    %make_dummy_tile_info create a info structure containing tile info
    %
    % use MATLAB georefpostings to return an info structure of type
    % map.rasterref.GeographicPostingsReference
    % Note the specific syntax here (set via spacing, not raster
    % size) is required to get a binary equivalent to what the
    % copernicus tiles return via readgeoraster(). E.g., this
    % method should avoid any small floating point round errors.
    % returns the single ref object.
        tlat = str2double(tile_id(2:3));
        tlon = str2double(tile_id(6:8));
        if tile_id(1) == 'S'
            tlat = tlat * -1;
        end
        if tile_id(5) == 'W'
            tlon = tlon * -1;
        end

        latspacing = 1/1200;
        lonspacing = elevation_model.copernicus_lonspacing(tlat);
        latlim = [tlat + latspacing, tlat+1];
        lonlim = [tlon, tlon+1 - lonspacing];

        info = georefpostings( ...
            latlim, lonlim, latspacing, lonspacing, ...
            'ColumnsStartFrom', 'north');
    end

    function [tile_ul, tile_lr] = compute_tile_coverage(info)
    %compute_tile_coverage determine the UL and LR corners of the actual raster coverage
    %
    % the actual raster coverage is computed based on the
    % assumption that the raster is defined as "postings", which
    % means the lat/lon positions can be viewed as the center of a
    % grid cell that extends 0.5x the sampline spacing in all
    % directions. Thus, the tile coverage between neighboring
    % tiles will have no gaps, while the "Limits" in the geo info
    % structure will have a gap equal to one SampleSpacing.
    %
    % returns tile_ul, tile_lr as 2 element arrays.
    % UL = upper left = NW corner; LR = lower right = SE corner.
    %
    % This is not currently used, but I've left it implemented in
    % case it is needed again.
        tile_ul = [ ...
            info.LongitudeLimits(1) - 0.5 * info.SampleSpacingInLongitude, ...
            info.LatitudeLimits(2)  + 0.5 * info.SampleSpacingInLatitude];
        tile_lr = [ ...
            info.LongitudeLimits(2) + 0.5 * info.SampleSpacingInLongitude, ...
            info.LatitudeLimits(1)  - 0.5 * info.SampleSpacingInLatitude];
    end

    function P = compute_tile_polyshape(info)
    %compute_tile_polyshape create a polyshape object based on a geo info object
    %
    % the polyshape includes the 0.5x sample spacing around the
    % postings (see compute_tile_coverage description)
    % returns one polyshape object
        lon = zeros(4,1);
        lat = zeros(4,1);
        % MATLAB polyshape prefers CW orientation of the corner points
        % UL -> UR -> LR -> LL, or, NW -> NE -> SE -> SW.
        lon(1) = info.LongitudeLimits(1) - 0.5 * info.SampleSpacingInLongitude;
        lon(2) = info.LongitudeLimits(2) + 0.5 * info.SampleSpacingInLongitude;
        lon(3) = info.LongitudeLimits(2) + 0.5 * info.SampleSpacingInLongitude;
        lon(4) = info.LongitudeLimits(1) - 0.5 * info.SampleSpacingInLongitude;

        lat(1) = info.LatitudeLimits(2) + 0.5 * info.SampleSpacingInLatitude;
        lat(2) = info.LatitudeLimits(2) + 0.5 * info.SampleSpacingInLatitude;
        lat(3) = info.LatitudeLimits(1) - 0.5 * info.SampleSpacingInLatitude;
        lat(4) = info.LatitudeLimits(1) - 0.5 * info.SampleSpacingInLatitude;

        P = polyshape(lon, lat);
    end

    function [polygon_ul, polygon_lr] = compute_polygon_coverage(P)
    %compute_polygon_coverage compute the UL, LR bounding corner points for the input polyshape
        [lon_range, lat_range] = boundingbox(P);
        polygon_ul = [lon_range(1), lat_range(2)];
        polygon_lr = [lon_range(2), lat_range(1)];
    end

    function [eastP, westP] = poly_antimeridian_unfold(P)
    %poly_antimeridian_unfold "twin" a polyshape into an East and West polyshape
    %
    % The input polyshape is twinned into an East polyshape, by
    % adding 360 degrees to negative longitudes, so the polyshape
    % will be continuous in the range [0, 360] degrees, and will
    % be valid for overlaps with E179 tiles;
    % and a West polyshape, by subtracting 360 from positive
    % longitudes, so the polyshape will be continuous in the range
    % [-360, 0] degrees, and will be valid for overlaps with W180
    % tiles.
    %
    % the input polyshape must contain the antimeridian, otherwise
    % bad values will be computed.
        lon = P.Vertices(:,1);
        lat = P.Vertices(:,2);
        unfolded_east_lon = lon;
        unfolded_east_lon(lon < 0) = lon(lon < 0) + 360.0;
        unfolded_west_lon = lon;
        unfolded_west_lon(lon > 0) = lon(lon > 0) - 360.0;
        eastP = polyshape(unfolded_east_lon, lat);
        westP = polyshape(unfolded_west_lon, lat);
    end


    function stat = contains_antimeridian(P)
    % "Walk" along the outer boundary points of a polygon.
    %
    % Polygons that cross neither the antimeridian nor a pole will
    % have no sudden discontinuities in longitude values.
    %
    % Polygons that cross the antimeridian but not a pole will have
    % two discontinuities.
    %
    % Polygons that cross a single pole will have one
    % discontinuity.
    %
    % Since this framework is assumed to work only on TIRS scene
    % footprints, we make the assumption that there are no polygons
    % that contain the pole. If that is the case, then we throw an
    % exception as that case will fail in unexpected ways.
    %
    % therefore, the return 'stat' is true (contains antimeridian)
    % or false (does not contain antimeridian)
    %
    % Note: this is adapted from the similar function in
    % PREFIRE_AUX_SAT, written by K.

        % Note: this method won't capture the half-pixel on the
        % other side of the antimeridian from the tiles on the East
        % side. (e.g., the W180 tile, which borders the
        % antimeridian to the east, includes 0.5*lonspacing on the
        % west side of the antimerdian.) I think, because the
        % overlapping is then applied the the raster postings, this
        % shouldn't cause any real issue.

        [boundary_lon, ~] = boundary(P);
        lon_abs_diffs = abs(diff(boundary_lon));
        discts = lon_abs_diffs > 180.0;
        num_crosses = sum(discts);
        if num_crosses == 1
            error = MException("elevation_model:GeometryError", ...
                "polygon appears to contain a pole");
            throw(error)
        end
        stat = num_crosses == 2;
    end

    function [P1, P2] = poly_antimeridian_split(P)
    % Split polygons that cross the antimeridian into two polygons
    % that have the -180 longitude line as one edge.
    % This function assumes the polygon does cross the antimeridian.
    %
    % Note: this is adapted from the similar function in
    % PREFIRE_AUX_SAT, written by K.
    %
    % This is not currently used, as the 'twinning' or 'unfolding'
    % method was implemented instead of this.
    % (see poly_antimeridian_unfold)
    %
    % Originally, I believed this splitting method
    % would be intractable because the actual border between the
    % tiles is 0.5x pixelspacing to the west of the antimeridian,
    % so splitting at exactly 180 degrees doesn't split the scene
    % polygon to match the tile edges.
    % However, since the overlapping calculation is applied to the
    % raster postings (which do effectively split right at 180.0),
    % that offset might not matter. Leaving this method implemented
    % in case we need it later, if the twinning doesn't seem to
    % work for some corner cases I had not noticed.
    % For this splitting to work, I think it would need to split at
    % +179.99999 longitude (for the East tile) so that it doesn't
    % try to include the actual antimeridian; the copernicus tiles
    % contain the antimerdian in the W180 tiles at -180.0000 deg.

        [bdry_lon, bdry_lat] = boundary(P);

        % compute differences in longitude on boundary:
        % keep track of East to West, and West to each jumps (from
        % the perspective of "walking" the boundary.)
        londiff = diff(bdry_lon);
        split_ixs = find(abs(londiff) > 180);
        E2Wfirst = londiff(split_ixs(1)) < -180;

        % Finally, generate list of points for each of the split polygons.
        % One edge of each polygon is created by drawing boundary directly
        % along the -180 / 180 longitude line.
        % The latitudes of these edge points are calculated as the mean
        % latitude of the two points on either side of the split.
        % TBD: this should be (probably) the linear interpolation,
        % rather than the average. matters if the polygon is 'large'.
        % FIXME
        split_lat_1 = ( bdry_lat(split_ixs(1)) + bdry_lat(split_ixs(1)+1) ) / 2;
        split_lat_2 = ( bdry_lat(split_ixs(2)) + bdry_lat(split_ixs(2)+1) ) / 2;

        % note that the copernicus tiles start at the NW corner,
        % and do not include postings at the E edge. So, the
        % inserted split values need to be just inside the
        % antimeridian, by the longitude increment. Here,
        % we've effectively hardcoded in the grid spacing, which is
        % not great, but it would be messy to get the georaster info
        % structure here.

        %
        % note on lon increment: the changing lon spacing makes
        % this tricky. the problem here is that the W180 tiles have
        % converage by 0.5*lonspacing into the E hemisphere. If
        % this polygon cross a latitude where the lonspacing
        % changes, then we need to choose the coarser lonspacing in
        % order to ensure the split poly doesn't just cross the
        % edge of the tile on the E side of the Antimeridian (the
        % W180 tile.)
        min_lat = min(bdry_lat);
        tlat = elevation_model.lat_to_copernicus_tilelat(min_lat);
        lonspacing1 = elevation_model.copernicus_lonspacing(tlat);

        max_lat = max(bdry_lat);
        tlat = elevation_model.lat_to_copernicus_tilelat(max_lat);
        lonspacing2 = elevation_model.copernicus_lonspacing(tlat);

        lonspacing = max([lonspacing1 lonspacing2]);

        if E2Wfirst
            lon1 = [bdry_lon(1:split_ixs(1)); ...
                    180.0 - lonspacing; 180.0 - lonspacing; ...
                    bdry_lon(split_ixs(2)+1:end)];
            lat1 = [bdry_lat(1:split_ixs(1)); ...
                    split_lat_1; split_lat_2; ...
                    bdry_lat(split_ixs(2)+1:end)];

            lon2 = [-180.0; ...
                    bdry_lon(split_ixs(1)+1:split_ixs(2)); ...
                    -180.0];
            lat2 = [split_lat_1; ...
                    bdry_lat(split_ixs(1)+1:split_ixs(2)); ...
                    split_lat_2];
        else
            lon1 = [bdry_lon(1:split_ixs(1)); ...
                    -180.0; -180.0; ...
                    bdry_lon(split_ixs(2)+1:end)];
            lat1 = [bdry_lat(1:split_ixs(1)); ...
                    split_lat_1; split_lat_2; ...
                    bdry_lat(split_ixs(2)+1:end)];

            lon2 = [180.0 - lonspacing; ...
                    bdry_lon(split_ixs(1)+1:split_ixs(2)); ...
                    180.0 - lonspacing];
            lat2 = [split_lat_1; ...
                    bdry_lat(split_ixs(1)+1:split_ixs(2)); ...
                    split_lat_2];
        end

        P1 = polyshape(lon1, lat1);
        P2 = polyshape(lon2, lat2);

    end

    function ll_grid = make_ll_grid(tile_info)
    %make_ll_grid compute the per-pixel longitude and latitude arrays for the tile's raster postings.
    %
    % the created array is shaped (RasterSize(1), RasterSize(2), 2),
    % where (:,:,1) are the longitudes and (:,:,2) are the
    % latitudes. The positions are the raster postings, which I
    % assume to be the centers of grid cells with SampleSpacing size.
        ll_grid = zeros([tile_info.RasterSize(1), tile_info.RasterSize(2), 2]);
        lon_vec = tile_info.SampleSpacingInLongitude * ...
                  (0:tile_info.RasterSize(2)-1) + ...
                  tile_info.LongitudeLimits(1);
        lat_vec = tile_info.SampleSpacingInLatitude * ...
                  (tile_info.RasterSize(1)-1:-1:0)' + ...
                  tile_info.LatitudeLimits(1);
        ll_grid(:,:,1) = repmat(lon_vec, [tile_info.RasterSize(1), 1]);
        ll_grid(:,:,2) = repmat(lat_vec, [1, tile_info.RasterSize(2)]);
    end

    function [v1, v2] = extract_pixel_group(...
        tile1, tile_ll, tile_info, P, tile2)
    %extract_pixel_group core low level function to extract DEM pixels
    %
    % inputs are:
    % tile1: raster array for extraction of pixel values, with
    %     shape RasterSize, (i,j) 
    % tile_ll the lon-lat per raster pixel, shaped (i,j,2) to match
    %     the tile data in tile1.
    % tile_info: the geo info structure, created by readgeoraster.
    % P: the polyshape describing the region to extract
    % tile2: an optional second raster array for extraction. Must
    %     match tile1's shape.
    %
    % tile1 and tile2 can be arrays of any data type or contents,
    % as long as the have the correct shapes and match tile_ll.
    % By implementing tile2 as optional we cover two use cases:
    % when two arrays are extracted at the same time with P. This
    % is the typical case (extract altitude and mask, or lon and
    % lat), but then for an ocean tile we only need to do one
    % lookup (to know the overlap size) as the mask is known.

        % tile 2 is optional
        extract_tile2 = nargin == 5;

        % extract bounding subset around the scene poly
        % [j1,i1] is the lon_j, lat_i for Upper Left corner;
        % [j2,i2] is for the Lower Right.
        [polygon_ul, polygon_lr] = elevation_model.compute_polygon_coverage(P);
        [j1,i1] = elevation_model.ll_to_ij(tile_info, polygon_ul(1), polygon_ul(2));
        [j2,i2] = elevation_model.ll_to_ij(tile_info, polygon_lr(1), polygon_lr(2));

        % expand by one tile pixel to give a buffer
        if i1 > 1
            i1 = i1 - 1;
        end
        if j1 > 1
            j1 = j1 - 1;
        end
        if i2 < tile_info.RasterSize(1)
            i2 = i2 + 1;
        end
        if j2 < tile_info.RasterSize(2)
            j2 = j2 + 1;
        end

        %disp(sprintf('extract_h_group, corner pix: %5d %5d %5d %5d',i1,j1,i2,j2));
        tile1_sub = tile1(i1:i2,j1:j2);
        lon_sub = tile_ll(i1:i2,j1:j2,1);
        lat_sub = tile_ll(i1:i2,j1:j2,2);
        nvals = numel(lat_sub);

        lon_sub_1d = reshape(lon_sub, [nvals,1]);
        lat_sub_1d = reshape(lat_sub, [nvals,1]);
        %disp(sprintf('extract_h_group, lon/lat min: %15.10f %15.10f', ...
        %             min(lon_sub_1d), min(lat_sub_1d)));
        %disp(sprintf('extract_h_group, lon/lat max: %15.10f %15.10f', ...
        %             max(lon_sub_1d), max(lat_sub_1d)));

        msk_1d = isinterior(P, lon_sub_1d, lat_sub_1d);
        msk = reshape(msk_1d, size(tile1_sub));

        v1 = tile1_sub(msk);
        if extract_tile2
            tile2_sub = tile2(i1:i2,j1:j2);
            v2 = tile2_sub(msk);
        end

    end


    function dlon = copernicus_lonspacing(tile_lat)
    % copernicus_lonspacing compute the longitude spacing give the tile latitude
    %
    % here, tile latitude is the integer latitude number that is
    % also part of the tile's filename. this will not work properly
    % for arbitrary longitude values, for values near tile borders.
        if tile_lat >= 85
            dlon = 1/120;
        elseif tile_lat >= 80
            dlon = 1/240;
        elseif tile_lat >= 70
            dlon = 1/400;
        elseif tile_lat >= 60
            dlon = 1/600;
        elseif tile_lat >= 50
            dlon = 1/800;
        elseif tile_lat >= -50
            dlon = 1/1200;
        elseif tile_lat >= -60
            dlon = 1/800;
        elseif tile_lat >= -70
            dlon = 1/600;
        elseif tile_lat >= -80
            dlon = 1/400;
        elseif tile_lat >= -85
            dlon = 1/240;
        else
            dlon = 1/120;
        end
    end

    function [tlat] = lat_to_copernicus_tilelat(lat)
    %lat_to_copernicus_tilelat convert an arbitrary fractional latitude value to integer tile latitude
    %
    % the returned tile lat is an integer values, which is used
    % in the filename construction for the tile. The calculation
    % ensures the tile will contain the requested arbitrary
    % fractional lat position.
        dlat = 1/1200;
        tlat = floor(lat - 0.5*dlat);
    end

    function [tlon, tlat] = ll_to_copernicus_tilell(lon, lat, offset)
    % ll_to_copernics_tilell return integer tile lon/lat for the input fractional lon/lat.
    %
    % the returned tile lon, lat, are integer values, that are used
    % in the filename construction for the tile. The calculation
    % ensures the tile will contain the requested arbitrary
    % fractional lon/lat position.

        % here, convert lat to tile lat first (this has a fixed
        % spacing of 1/1200), as the lon spacing depends on the
        % tile latitude.
        % method here, puts the 'boundary' between tiles at half of the
        % spacing increment; this matches how
        % compute_tile_coverage() works.
        if nargin == 3
            lonspacing_scale = offset;
        else
            lonspacing_scale = 0.5;
        end
        tlat = elevation_model.lat_to_copernicus_tilelat(lat);
        dlon = elevation_model.copernicus_lonspacing(tlat);
        tlon = floor(lon + lonspacing_scale*dlon);
        if tlon >= 180
            tlon = tlon - 360;
        end
    end

    function [tlon, tlat] = ll_to_tilell(tile_info, lon, lat)
    %ll_to_tilell convert lon/lat values to tile lon/lat
    %
    % convert arbitrary lat lon to tile lat lon - this is done
    % to ensure the tile contains the input lat lon.
    % in this case, we use the tile info directly (faster than
    % calling ll_to_copernicus_tilell, I think.)
    %
    % This is currently unused, but could be used in cases where
    % the tile_info is available. Most of the time we are using the
    % lat/lon to derive the tile lat/lon, though, in order to
    % figure out what tile to load. At that point we don't have the
    % tile info structure to use.
        dlon = tile_info.SampleSpacingInLongitude;
        dlat = tile_info.SampleSpacingInLatitude;
        tlon = floor(lon + 0.5*dlon);
        tlat = floor(lat - 0.5*dlat);
    end

    function tile_id = tilell_to_id(tlon, tlat)
    %tilell_to_id convert tile lon and lat to string tile_id
    %
    % the input tile lon and lat are convert to a string of the
    % form: (N|S)ii_(E|W)jjj where ii are the latitude digits and
    % jjj are the longitude digits.
        latstr = sprintf('%02d', abs(floor(tlat)));
        if tlat < 0
            lathalf = 'S';
        else
            lathalf = 'N';
        end
        lonstr = sprintf('%03d', abs(floor(tlon)));
        if tlon < 0
            lonhalf = 'W';
        else
            lonhalf = 'E';
        end
        tile_id = [lathalf latstr '_' lonhalf lonstr];

    end

    function [lon_j, lat_i] = ll_to_ij(tile_info, lon, lat)
    %ll_to_ij convert fractional lon and lat values to row/column indices
    %
    % the conversion is performed based on the parameters in
    % tile_info. Note that the i, j index values are forced to be
    % in valid range for the RasterSize described in the
    % tile_info. This means the returned lon_j, lat_i should never
    % be out of range index values. This also implies that an input
    % location that was actually outside the tile bounds would be
    % forced to the nearest tile edge point.

        % assumes even size grid cells in lon,lat ('equal angle')
        % note that the (1,1) element is the northwest corner, and
        % that lat is the first index (row number) in the 2D
        % array. Hence, the position in tile T will be T(lat_i,lon_j)
        dlon = lon - tile_info.LongitudeLimits(1);
        dlat = tile_info.LatitudeLimits(2) - lat;
        lon_j = round(dlon / tile_info.SampleSpacingInLongitude);
        lat_i = round(dlat / tile_info.SampleSpacingInLatitude);
        if lat_i <= 0
            lat_i = 1;
        end
        if lon_j <= 0
            lon_j = 1;
        end
        if lat_i > tile_info.RasterSize(1)
            lat_i = tile_info.RasterSize(1);
        end
        if lon_j > tile_info.RasterSize(2)
            lon_j = tile_info.RasterSize(2);
        end
    end

    function h = nearest_neighbor_h(tile, tile_info, lon, lat)
    %nearest_neighbor_h find tile data at posting nearest to input lon,lat position
    %
    % using tile_info to convert the lon lat to tile index i,j,
    % then return the tile value at that point.
        [lon_j, lat_i] = elevation_model.ll_to_ij(tile_info, lon, lat);
        h = tile(lat_i, lon_j);
    end
        
end

end
