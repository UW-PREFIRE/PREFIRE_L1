% redo the steps in read_REFIRE_payload_L0 between 'srd' array and
% 'science frame' arrays.

%%%
% forward path: payload SRD to Sci

% create dummy SRD frames where the value is equal to the position
% index in the raw packet frame (1 - 512)

addpath('../functions');

nframe = 3;
srd_counts_f = ones(512,nframe);
for n=1:nframe
    srd_counts_f(:,n) = (1:512)';
end

% reorder from (512,N) to (64,8,N) and the polarity flip is done in
% the same single loop.
srd_tmp_f = zeros(64,8,nframe);
for j=1:8
    for k=1:64
        srd_tmp_f(k,j,:) = ((-1)^k) * srd_counts_f(8*(k-1)+j,:);
    end
end

% apply ROIC map.
sci_counts_f = ROIC_map_io(64,8,nframe,srd_tmp_f);

%%%
% reverse path: Sci to Payload SRD
sci_counts_r = sci_counts_f * 1; % copy in MATLAB?
for j=1:4
    for k=1:64
        sci_counts_r(k,j,:) = (-1)^k * sci_counts(k,j,:);
    end
end
for j=5:8
    for k=1:64
        sci_counts_r(k,j,:) = (-1)^(k) * sci_counts(k,j,:);
    end
end

srd_tmp_r = ROIC_map_oi([8,nframe],64,sci_counts_r);
srd_counts_r = zeros(512,nframe);
for j=1:8
    for k=1:64
        srd_counts_r(8*(k-1)+j,:) = srd_tmp_r(k,j,:);
    end
end

save('ROIC_polarity_test.mat', ...
     'srd_counts_f', 'srd_tmp_f', 'sci_counts_f', ...
     'srd_counts_r', 'srd_tmp_r')


