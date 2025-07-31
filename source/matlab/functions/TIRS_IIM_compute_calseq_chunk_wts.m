function [wt] = TIRS_IIM_compute_chunk_wts(frame_ctimes, calseq_ctimes)
%
% [wt] = TIRS_IIM_compute_chunk_wts(frame_ctimes, calseq_ctimes)
%
% Compute weighting array, shaped (nframe, 4), across a single
% calseq data chunk.
% Frame_ctimes is assumed to be shaped (1,nframe+1), which is the
% minimum frame set to enclose the calseq_ctimes; meaning,
% calseq_ctimes(1) is between frame_ctimes(1) and frame_ctimes(2),
% and calseq_ctimes(end) is between frame_ctimes(end-1) and
% frame_ctimes(end).
%
% frame types
% 1 = earth
% 2 = space
% 3 = internal cal target
% 4 = instrument housing

% calseq_ctimes:
%
%     earth1 view
% 0 earth1 -> instr
%     instrument view
% 1 instr -> space
%     space view
% 2 space -> instr
%     instrument view
% 3 instr -> calib
%     calib view
% 4 calib -> instr
%     instrument view
% 5 instr -> earth2
%     earth2 view

nframes = length(frame_ctimes) - 1;
ntypes = 4;

wt = zeros([nframes, ntypes]);
frame_durations = diff(frame_ctimes);
calseq_durations = diff(calseq_ctimes);

% assumptions:
% start of calseq should always occur in frame 1.
% first calseq segment (the instrument view) is shorter than 1 frame time.

% first frame.
wt(1,1) = (calseq_ctimes(1) - frame_ctimes(1)) / frame_durations(1);
if calseq_ctimes(2) < frame_ctimes(2)
    wt(1,4) = calseq_durations(1) / frame_durations(1);
    wt(1,2) = (frame_ctimes(2) - calseq_ctimes(2)) / frame_durations(1);
else
    wt(1,4) = (frame_ctimes(2) - calseq_ctimes(1)) / frame_durations(1);
end

% second frame.
% assumptions: calseq segment (space view) is long (> 2 frames?)
if calseq_ctimes(2) > frame_ctimes(2)
    wt(2,4) = (calseq_ctimes(2) - frame_ctimes(2)) / frame_durations(2);
    wt(2,2) = (frame_ctimes(3) - calseq_ctimes(2)) / frame_durations(2);
else
    wt(2,2) = 1.0;
end

% whole space view frames
i1 = single_searchsorted(frame_ctimes, calseq_ctimes(3)) - 1;
wt(3:i1-1,2) = 1.0;

% mixed space / instrument, possible mix to calib, first frame
wt(i1,2) = (calseq_ctimes(3) - frame_ctimes(i1)) / frame_durations(i1);
if calseq_ctimes(4) < frame_ctimes(i1+1)
    wt(i1,4) = calseq_durations(3) / frame_durations(i1);
    wt(i1,3) = (frame_ctimes(i1+1) - calseq_ctimes(4)) / frame_durations(i1);
else
    wt(i1,4) = (frame_ctimes(i1+1) - calseq_ctimes(3)) / frame_durations(i1);
end

% possible mixed instr / calib , second frame
if calseq_ctimes(4) > frame_ctimes(i1+1)
    wt(i1+1,4) = (calseq_ctimes(4) - frame_ctimes(i1+1)) / frame_durations(i1+1);
    wt(i1+1,3) = (frame_ctimes(i1+2) - calseq_ctimes(4)) / frame_durations(i1+1);
else
    wt(i1+1,3) = 1.0;
end

% whole calib view frames
% assumptions: calseq segment (calib view) is long (> 2 frames?)
i2 = single_searchsorted(frame_ctimes, calseq_ctimes(5)) - 1;
wt(i1+2:i2,3) = 1.0;

% mixed calib / instrument, first frame
wt(i2,3) = (calseq_ctimes(5) - frame_ctimes(i2)) / frame_durations(i2);
if calseq_ctimes(6) < frame_ctimes(i2+1)
    wt(i2,4) = calseq_durations(5) / frame_durations(i2);
    wt(i2,1) = (frame_ctimes(i2+1) - calseq_ctimes(6)) / frame_durations(i2);
else
    wt(i2,4) = (frame_ctimes(i2+1) - calseq_ctimes(5)) / frame_durations(i2);
end

% possible mixed instr /earth 2, second frame
if calseq_ctimes(6) > frame_ctimes(i2+1)
    wt(i2+1,4) = (calseq_ctimes(6) - frame_ctimes(i2+1)) / frame_durations(i2+1);
    wt(i2+1,1) = (frame_ctimes(i2+2) - calseq_ctimes(6)) / frame_durations(i2+1);
end
% nothing to do here for else condition.
% this would only arise if there where more frames to process,
% and we assume the caller gave us only this limited chunk.

end


function [idx] = single_searchsorted(x, s)
% helper to mimic the operation of np.searchsorted.
% input x is assumed to be a sorted 1D vector, and s is a scalar.
% returns the index value, where s would be inserted into x,
% keeping x sorted.
[~,sort_idx] = sort([x, s]);
num = length(x);
idx = find(sort_idx == (num+1));
end
