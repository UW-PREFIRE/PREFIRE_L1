function calseqs = TIRS_IIM_create_calseq_chunks( ...
    calseqs, ch_irad, ctime, BRF, T_cal, T_instr);
%
% [calseqs] = TIRS_IIM_create_calseq_chunks(calseqs, ch_irad, BRF, T_cal, T_instr)
%
% Extract frame chunks from ch_irad (assumed shape (64,8,nframe))
% calseqs structure. Since the chunks may vary by +/- 1 frame in
% size, the chunks need to be stored in a cell array.
%
% the chunks of irad data will contain the calseq radiances,
% linearly mixed as appropriate, according to the detailed cal seq
% timing.
% each chunk is shaped (64,8,nframes_in_calseq)
%
% Inputs:
% calseqs: structure from TIRS_IIM_prep_calseqs()
% ch_irad: input channel integrated radiance array, shaped
%     (64,8,nframe). Normally this is from simulated data.
% ctime: array of ctime associated with the ch_irad, shaped (nframe)
% BRF: background response function data structure, loaded from stored
%     MATLAB file (this file is loaded within
%     TIRS_inverse_instrument_model())
% T_cal: time series of internal cal target temperature. Should be
%     a 1D vector, with nframe elements.
% T_instr: time series of internal instrument temperatures. Should
%     be a 1D vector, with nframe elements (same as T_cal).
%
% outputs:
% calseq structure, now including irad_chunks, a cell array with
%     num_calseq elements.
%     each irad_chunk has (64,8,nframe) elements, where nframe is
%     the number of frames in the single chunk.
%

[nwave, nscene, nframe_in] = size(ch_irad);
calseqs.irad_chunks = cell([1,calseqs.num]);

for n = 1:calseqs.num

    % First, load chunks into a (64,8,nframe,4) array, with the
    % additional axis containing the four different types:
    % 1 = Earth (Nadir), 2 = Space
    % 3 = internal calibration target, 4 = instrument housing

    a = calseqs.start_frames(n);
    b = calseqs.end_frames(n);
    nframe = b - a + 1;
    irad_chunk_t = zeros([nwave,nscene,nframe,4]);
    % slot 1 = earth view, copy all frames for the chunk
    irad_chunk_t(:,:,:,1) = ch_irad(:,:,a:b);
    % slot 2 = space view, constant, so fill all frames in the chunk
    irad_chunk_t(:,:,:,2) = ...
        repmat(BRF.scene_irad_s, [1,1,nframe]);
    % for other two views, loop through frames, as the temperatures
    % may be changing.
    for i = a:b
        % slot 3 = cal target view
        dT = T_cal(i) - BRF.T_cal(1);
        fraci = mod(dT, 1);
        i_cT = fix(dT)+2;
        irad_chunk_t(:,:,i-a+1,3) = ...
            (1-fraci) * BRF.scene_irad_c(:,:,i_cT-1) + ...
               fraci  * BRF.scene_irad_c(:,:,i_cT);
        % slot 4 = instrument internal view
        dT = T_instr(i) - BRF.T_cal(1);
        fraci = mod(dT, 1);
        i_cT = fix(dT)+2;
        irad_chunk_t(:,:,i-a+1,4) = ...
            (1-fraci) * BRF.scene_irad_c(:,:,i_cT-1) + ...
               fraci  * BRF.scene_irad_c(:,:,i_cT);
    end

    % Get weights from helper function
    idx = find( (ctime > calseqs.start_ctimes(n)) & ...
                (ctime <= calseqs.end_ctimes(n)) );
    frame_ctime_n = ctime(idx(1)-1:idx(end)+1);
    wt = TIRS_IIM_compute_calseq_chunk_wts(frame_ctime_n, calseqs.event_ctimes(:,n));

    % apply weights
    irad_chunk = zeros([nwave,nscene,nframe-1]);
    for f = 1:nframe-1
        for s = 1:4
            irad_chunk(:,:,f) = irad_chunk(:,:,f) + ...
                irad_chunk_t(:,:,f,s) * wt(f,s);
        end
    end

    calseqs.irad_chunks{n} = irad_chunk;
            
end
 
end
