function [signew] = ROIC_map_oi(num_images,spectral_channels,sig);
num_xtrack = num_images(1);
num_intrack = num_images(2);
% incoming signal array is ordered from science processing
% outgoing signal array is raw data for TIRS focal plane simulation

% incoming signal array is 8x64x(intrack) with the 8 dimension mapping into
% rows A-H and the 64 dimension mapping into channels (called pixels in Kenyon's documentation) 0-63,  the ROIC
% mapping is 1:1 for Rows A,B and G,H, but cut halfsies for Rows C-F

         
%read order    A0,   A1,   B0,   B1,   C0,   C1,   D0,  D1 (ROIC/CORE)
%               1     2     3     4     5     6     7     8   
%science order AS,AL,BS,BL,CS,CL,DS,DL,ES,EL,FS,FL,GS,GL,HS,HL

signew = zeros([spectral_channels,num_xtrack,num_intrack]);

%map order      <>  science order
signew( 1:32,1,:) =      sig( 1:32,4,:);   %A0a-DS
signew(33:64,1,:) =      sig( 1:32,3,:);   %A0b-CS
signew( 1:32,2,:) = flip(sig( 1:32,6,:),1);%A1a-FS r
signew(33:64,2,:) = flip(sig( 1:32,5,:),1);%A1b-ES r
signew( 1:32,3,:) =      sig(33:64,1,:);   %B0a-AL
signew(33:64,3,:) =      sig(33:64,2,:);   %B0b-BL
signew( 1:32,4,:) =      sig( 1:32,2,:);   %B1a-BS
signew(33:64,4,:) =      sig( 1:32,1,:);   %B1b-AS
signew( 1:32,5,:) = flip(sig(33:64,5,:),1);%C0a-EL r
signew(33:64,5,:) = flip(sig(33:64,6,:),1);%C0b-FL r
signew( 1:32,6,:) =      sig(33:64,3,:);   %C1a-CL
signew(33:64,6,:) =      sig(33:64,4,:);   %C1b-DL
signew( 1:32,7,:) = flip(sig( 1:32,8,:),1);%D0a-HS r
signew(33:64,7,:) = flip(sig( 1:32,7,:),1);%D0b-GS r
signew( 1:32,8,:) = flip(sig(33:64,7,:),1);%D1a-GL r
signew(33:64,8,:) = flip(sig(33:64,8,:),1);%D1b-HL r

end
