function [srd_sci] = ROIC_map_io(n_spectral, n_xtrack, n_atrack, srd);

% Incoming 'srd' array is raw (but reordered into a 3-D array, with the
%  alternating channel DN polarity) data from a TIRS focal plane.
%
% Outgoing 'srd_sci' array is ordered for science processing.
%
% 8x64x(along-track) with the xtrack dimension mapping into rows A-H. The ROIC
%  mapping is 1:1 for Rows A,B and G,H, but cut halfsies for Rows C-F.
%
% read order    A0,   A1,   B0,   B1,   C0,   C1,   D0,   D1  (ROIC:CORE)
%               1     2     3     4     5     6     7     8   
% science order AS,AL,BS,BL,CS,CL,DS,DL,ES,EL,FS,FL,GS,GL,HS,HL

srd_sci = zeros([n_spectral,n_xtrack,n_atrack]);

% science order      <>  detector order
srd_sci( 1:32,4,:) =      srd( 1:32,1,:);     %DS-A0a
srd_sci( 1:32,3,:) =      srd(33:64,1,:);     %CS-A0b
srd_sci( 1:32,6,:) = flip(srd( 1:32,2,:), 1); %FS-A1a r
srd_sci( 1:32,5,:) = flip(srd(33:64,2,:), 1); %ES-A1b r
srd_sci(33:64,1,:) =      srd( 1:32,3,:);     %AL-B0a
srd_sci(33:64,2,:) =      srd(33:64,3,:);     %BL-B0b
srd_sci( 1:32,2,:) =      srd( 1:32,4,:);     %BS-B1a
srd_sci( 1:32,1,:) =      srd(33:64,4,:);     %AS-B1b
srd_sci(33:64,5,:) = flip(srd( 1:32,5,:), 1); %EL-C0a r
srd_sci(33:64,6,:) = flip(srd(33:64,5,:), 1); %FL-C0b r
srd_sci(33:64,3,:) =      srd( 1:32,6,:);     %CL-C1a
srd_sci(33:64,4,:) =      srd(33:64,6,:);     %DL-C1b
srd_sci( 1:32,8,:) = flip(srd( 1:32,7,:), 1); %HS-D0a r
srd_sci( 1:32,7,:) = flip(srd(33:64,7,:), 1); %GS-D0b r
srd_sci(33:64,7,:) = flip(srd( 1:32,8,:), 1); %GL-D1a r
srd_sci(33:64,8,:) = flip(srd(33:64,8,:), 1); %HL-D1b r

end
