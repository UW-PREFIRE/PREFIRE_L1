function [error_ca, calseq] = process_TIRS_calseqs(L0, L1, sig)
% Determine which calibration sequences are valid, and process those further.
% The full L0 structure and a partial L1 structure are input. The
% L1 needs only to have the instrument temperatures populated.
%
% OUTPUT fields:
%   calseq.n:  Number of valid calibration sequences
%   calseq.space_seqmean_ts:  mean ctime of space frames in each calseq [s]
%   calseq.space_seqmean_v:  mean (across space frames in each calseq) values
%   calseq.contrast_seqmean_ts:  mean ctime of each calseq centroid [s]
%   calseq.contrast_seqmean_v:  mean contrast (for each calseq) values
%   calseq.cal_seqmean_ts:  mean ctime of cal target frames
%   calseq.cal_seqmean_temp:  mean calibrator temperature for the calseq [K]
%   calseq.space_framecount:  number of space frames in each calseq
%   calseq.cal_framecount:  number of cal target frames in each calseq
error_ca = {'#NONE#', 0};  % Default

% Determine the available full blocks of calibration-target frames:
%  - an edge is indicated by > 25 seconds between frames
ic_all_block_edges = find(diff(L0.caltgt.ctime) > 25);
if_cal_ld_edges = ic_all_block_edges(1:end-1)+1;
if_cal_tr_edges = ic_all_block_edges(2:end);
n_caltgt_blocks = length(if_cal_ld_edges);

% Determine the available full blocks of space frames:
%  - an edge is indicated by > 25 seconds between frames
ic_all_block_edges = find(diff(L0.space.ctime) > 25);
if_spc_ld_edges = ic_all_block_edges(1:end-1)+1;
if_spc_tr_edges = ic_all_block_edges(2:end);
n_space_blocks = length(if_spc_ld_edges);

% Perform some sanity checks:
if (n_caltgt_blocks ~= n_space_blocks)
   error_ca = {['ERROR: Number of usable space frame blocks must ' ...
                 'equal number of usable caltgt frame blocks.'], 21};
   calseq = -999;  % Allow error code propagation (avoid "not assigned")
   return
end
if (n_space_blocks == 0)
   error_ca = {'ERROR: no usable calibration sequences found.', 20};
   calseq = -999;  % Allow error code propagation (avoid "not assigned")
   return
end

% look through detected calseq and remove any that are abnormally
% long (these are often cal target stares during contacts, for example)
max_valid = 30;
max_num_blocks = n_caltgt_blocks;
valid_blocks = logical(ones([max_num_blocks,1]));
for n = 1:n_space_blocks
    space_view_length = if_spc_tr_edges(n) - if_spc_ld_edges(n) + 1;
    calib_view_length = if_cal_tr_edges(n) - if_cal_ld_edges(n) + 1;
    if (space_view_length > max_valid) | (calib_view_length > max_valid)
        valid_blocks(n) = false;
    end
end

n_valid_blocks = sum(valid_blocks);
if_cal_ld_edges = if_cal_ld_edges(valid_blocks);
if_cal_tr_edges = if_cal_tr_edges(valid_blocks);
if_spc_ld_edges = if_spc_ld_edges(valid_blocks);
if_spc_tr_edges = if_spc_tr_edges(valid_blocks);

calseq.n = n_valid_blocks;  % Number of usable calibration sequences

% Calculate 'space' and 'contrast' fields for all valid calibration sequences:
calseq.contrast_seqmean_ts = zeros(1, calseq.n);
calseq.contrast_seqmean_v = zeros(calseq.n, L0.pdims.allspectral, ...
                                  L0.pdims.xtrack);
calseq.space_seqmean_ts = zeros(1, calseq.n);
calseq.space_seqmean_v = zeros(calseq.n, L0.pdims.allspectral, L0.pdims.xtrack);
calseq.space_framecount = zeros(1, calseq.n);
calseq.cal_seqmean_ts = zeros(1, calseq.n);
calseq.cal_seqmean_v = zeros(calseq.n, L0.pdims.allspectral, L0.pdims.xtrack);
calseq.cal_framecount = zeros(1, calseq.n);
calseq.cal_seqmean_temp = zeros(1, calseq.n);
calseq.cal_seqmean_TIRS_T = zeros(8, calseq.n);
% along the way, as a simple filter, compute std/mean for each sequence, at
% spectral index = 15. This is one channel that is good in all 16 scenes.

% either one (space or contrast) is > 1% remove that sequence. This
% should be robust unless either uncalibrated signal resides around 0.
good_seq = logical(zeros(1,calseq.n));
test_ch = 15;
for i_cseq=1:calseq.n
   ib = if_spc_ld_edges(i_cseq);
   ie = if_spc_tr_edges(i_cseq);
   calseq.space_seqmean_ts(i_cseq) = mean(L0.space.ctime(ib:ie));  % [s]
   calseq.space_seqmean_v(i_cseq,:,:) = mean(sig.space(:,:,ib:ie), 3);
   calseq.space_framecount(i_cseq) = ie-ib+1;

   % MATLAB's array shape enforcement is weird: "squeezing" both does not
   % work, you need to transpose the first one because it is 2D.
   space_rel_std = std(sig.space(test_ch,:,ib:ie),0,3)' ./ ...
       squeeze(abs(calseq.space_seqmean_v(i_cseq,test_ch,:)));

   ib = if_cal_ld_edges(i_cseq);
   ie = if_cal_tr_edges(i_cseq);

   % Calculate spatial-mean temperature of the internal calibration target
   % (temps 4,5,6 are different locations in the calibrator.)
   calseq.cal_seqmean_TIRS_T(:,i_cseq) = mean(L1.caltgt.T_val(:,ib:ie), 2);
   calseq.cal_seqmean_temp(i_cseq) = mean(L1.caltgt.T_val(4:6,ib:ie), 'all');
   calseq.cal_seqmean_ts(i_cseq) = mean(L0.caltgt.ctime(ib:ie));  % [s]
   calseq.cal_seqmean_v(i_cseq,:,:) = mean(sig.caltgt(:,:,ib:ie), 3);
   calseq.cal_framecount(i_cseq) = ie-ib+1;

   calseq.contrast_seqmean_ts(i_cseq) = 0.5*( ...
       calseq.cal_seqmean_ts(i_cseq) + ...
       calseq.space_seqmean_ts(i_cseq) );  % [s]

   calseq.contrast_seqmean_v(i_cseq,:,:) = calseq.cal_seqmean_v(i_cseq,:,:)- ...
                                           calseq.space_seqmean_v(i_cseq,:,:);
   caltgt_rel_std = std(sig.caltgt(test_ch,:,ib:ie),0,3)' ./ ...
       squeeze(abs(calseq.cal_seqmean_v(i_cseq,test_ch,:)));
   good_seq(i_cseq) = all((space_rel_std < 0.01) & (caltgt_rel_std < 0.01));
end
% remove bad seq if any were identifed.
if any(~good_seq)
    calseq.contrast_seqmean_ts = calseq.contrast_seqmean_ts(1,good_seq);
    calseq.contrast_seqmean_v = calseq.contrast_seqmean_v(good_seq,:,:);
    calseq.space_seqmean_ts = calseq.space_seqmean_ts(1,good_seq);
    calseq.space_seqmean_v = calseq.space_seqmean_v(good_seq,:,:);
    calseq.space_framecount = calseq.space_framecount(1,good_seq);
    calseq.cal_seqmean_TIRS_T = calseq.cal_seqmean_temp(:,good_seq);
    calseq.cal_seqmean_temp = calseq.cal_seqmean_temp(1,good_seq);
    calseq.cal_seqmean_ts = calseq.cal_seqmean_ts(1,good_seq);
    calseq.cal_seqmean_v = calseq.cal_seqmean_v(good_seq,:,:);
    calseq.cal_framecount = calseq.cal_framecount(1,good_seq);
    calseq.n = sum(good_seq);
end

end
