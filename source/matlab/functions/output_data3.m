function output_data3(fileout,num_images,cadence,orbit_time,spectral_channels,housekeeping,encoder,sig);
% Forwared model version 4 sends a full load of housekeeping from the
% periodic thermal model, this made chunk boundary math obsolete here,
% remove cadence related input and code
%num is the number of thermal models points in track
num_xtrack = num_images(1);
num_intrack = num_images(2);

binID = fopen(fileout,'w');% this erases any prior file of same name
checksum = uint32(0);
%checksum = 0;
%byte 0 length 1
%packet_type = dec2hex(0); %switch to 1 for housekeeping
%fi(v,s,w,f) returns a fixed-point object with value v, Signed property value s, word length w, and fraction length f. Fraction length can be greater than word length or negative
imagepacket = 0;
engineering = 1;
packet_type = fi(imagepacket,0,1,2); %hd(2)
%fwrite(binID,imagepacket,'integer*1');
%byte 1 length 1 (2)
APID_LSB = 1; %0xDC; hd(1)
%byte 2 length 1 (3)
sequence = 0; %5:0 sequence count MSB, 7:6 grouping flags hd(3)
%byte 3 length 1 (4) 
pingpong = 0;% switch to 255 and back to 0 every other packet, set in loop hd(4)
%byte 4-5 length 2 (6)

overhead = 41; % this covers encoder, plus header 
packet_length = 2*(num_xtrack*(spectral_channels+1)) + overhead;%2* for number of bytes for each 'count' 1065 and counting hd(5)

% writing one packet with header is 1041 bytes, with 1024 for image data

%byte 6-9 length 4 (10)
time = 0; % will be reset in loop hd(6)
%byte 10-11 length 2 (12)
pps = 0; %milliseconds from last PPS hd(7)
checksum = checksum + pps;
%byte 12-15 length 4 (16)
cmd_count = 0; %command executed count,  hd(8)
checksum = checksum + cmd_count;
%byte 16-19 length 4 (20) 
last_cmd = 11; %hd(9)
checksum = checksum + last_cmd;
%byte 20-23 length 4 (24) 
interrupt = 0; %hd(10)
rotations = 90; %hd(11)
%byte N length 2 (26)
checksum = checksum + interrupt; %hd(12)
header_fmt = '%d%d%d%d%d%d%d%d%d%d%d';
thermal_model_data();

% last chance check on the signal array, making sure no NaN leaked through.
if any(isnan(sig), 'all')
    disp('Warning, Some NaN in signal array')
end

housekeeping_counts = inverse_tempest_temperature(housekeeping);
this_chunk = ' %d';
this_line ='%f %d %s';
housekeeping_fmt = ' %d';
housekeeping_length = 16; %need to adjust with header
checksum_so_far = checksum;
for m=1:num_intrack; %used to be num_temps
    %need to rewrite digital header here
    checksum = checksum_so_far;
%    housebin = sprintf('%s',fi(housekeeping_counts(mm,:),0,16,0));
    housebin = sprintf('%s',housekeeping_counts(1:8,m)); % indices swapped vs. v3 and counted by m now
%    chunkbin = sprintf('%s',fi(round(sig(:,:,m)),0,16,0));
    chunkbin = sprintf('%s',round(sig(:,:,m)));
    fwrite(binID,imagepacket+8*APID_LSB,'uint8'); %write first two bytes to file hd(1), hd(2)
    %checksum = checksum + imagepacket; %APID_LSB is not included in API checksum
    sequence = mod(m,cadence); %will be counting number within a chunk 
    fwrite(binID,sequence,'uint8'); %hd(3)
    pingpong = 255*mod(m,2);
    %pp = sprintf('%s',fi(pingpong,0,8,0));
    fwrite(binID,pingpong,'uint8'); %hd(4)
    checksum = checksum + pingpong;
    fwrite(binID,1067,'uint16'); %this is written LSB first hd(5) %fwrite(binID,1065,'uint16'); changed also 6/27/2020
    fwrite(binID,1000*orbit_time(m),'uint32'); %hd(6) retain the 1ms digit %B. expanded time from 16 to 32 bit on 6/27/2020
%    checksum = checksum + orbit_time(m);
    fwrite(binID,pps,'uint8'); %wasting another byte? hd(7)
    fwrite(binID,cmd_count,'uint32'); %4 bytes to match tempest hd(8)
    fwrite(binID,last_cmd,'uint32'); %4 bytes to match tempest hd(9)
    fwrite(binID,rotations,'uint32'); %4 bytes to match tempest hd(10)
    fwrite(binID,interrupt,'uint16'); %hd(11)
    fwrite(binID,encoder(m), 'uint16');
    fwrite(binID,housekeeping_counts(1:8,m), 'uint16');
    fwrite(binID,sig(:,:,m), 'uint16');
    %% To calculate the checksum of an API frame:
    %% Add all bytes of the packet, except the start delimiter 0x7E and the length (the second and third bytes).
    % B. also excluded packet_type and sequence, starting after length
    %% Keep only the lowest 8 bits from the result.
    checksum = bitshift(checksum,24); %checksum is defined as uint32 so left shift to put lowest 8 bits in MSB
    checksum = bitshift(checksum,-24); %now shift back to LSB 
    %% Subtract this quantity from 0xFF.
    checksum = 255 - checksum;
    fwrite(binID, checksum, 'uint8'); %hd(12)
    %test
    %fprintf('image # %d %s\n', m, imagepacket);
    %fprintf('image # %d %s\n', m, sequence);
    %fprintf('image # %d %s\n', m, pingpong);
    %fprintf('image # %d %s\n', m, 1000*orbit_time(m));
    %fprintf('image # %d %s\n', m, pps);
    %fprintf('image # %d %s\n', m, cmd_count);
    %fprintf('image # %d %s\n', m, last_cmd);
    %fprintf('image # %d %s\n', m, rotations);
    %fprintf('image # %d %s\n', m, interrupt);
    %fprintf('image # %d %s\n', m, encoder(m));
    %fprintf('image # %d %s\n', m, housebin);
    %fprintf('image # %d %s\n', m, chunkbin);
    %fprintf('image # %d %s\n', m, checksum);
    %fprintf('image # %d %s%s%s%s%s%s%s%s%s%s%s%s%s\n', m, imagepacket, sequence, pingpong, 1000*orbit_time(m), pps, cmd_count, last_cmd, rotations,interrupt, encoder(m), housebin, chunkbin, checksum);
end
fclose(binID);
% desire housekeeping packet every ~128 science or image packets 
% currently thermal model 255 seconds or every 360 packets
% image packet
%% header
%% encoder position ( 0,26,46,66 )*16384 / 360 + offset (0 is safe, 20 is nadir, 40 is internal target, 60 is space)
%% timestamp (use either 16 or 17 frames + overhead, 1 frame is 1/24 sec)
%% ROIC0,C0,P0 (16bits) (core 0)
%% ROIC0,C0,P1 (16bits)
%% ...
%% ROIC0,C1,P0 (16bits) (core 1)
%% ROIC0,C1,P1 (16bits)
%% ...
%% ROIC1,C0,P0 (16bits) (core 0)
%% ROIC1,C0,P1 (16bits)
%% ...
%% continue with 4 ROICS, each with 2 cores of 128 pixels
% housekeeping packet
%% header
%% 3.3V IC   (12bits)
%% VCORECDH
%% OBT1CDH
%% OBT2CDH
%% ROIC RTBA
%% ROIC RTBB
%% ROIC RTBC
%% ROIC RTBD
%% THERM0
%% THERM1
%% THERM2
%% THERM3
%% THERM4
%% THERM5
%% THERM6
%% THERM7

end
