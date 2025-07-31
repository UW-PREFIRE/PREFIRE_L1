function [obs_bitflags] = create_observation_bitflags( ...
    granule_beg_ts, granule_end_ts, obs_ctime, calseq_ctime, ...
    orbphase, orbphase_qinfo, calibrator_T, scipkt_ctime, ...
    ti, TIRS_ID, has_bus_bitflags, Cal_Artifacts)

% Create 1D array of bitflags for observation quality, depending on
% orbitphase ranges (denoting eclipse features), calibrator temp,
% and the time of valid cal sequences.
%
% Inputs:
% granule_beg_ts: ctime of granule beginning [s]
% granule_end_ts: ctime of granule ending [s]
% obs_ctime: 1D array of observation ctimes, with a number of
%     elements equal to the number of along-track data frames [s]
% calseq_ctime: 1D array of cal sequence ctimes [s]
% orbphase: 1D array of orbitphase values (valid at obs_ctime values) [deg]
% orbphase_qinfo: struct with orbitphase flagging information.
%     returned from read_PREFIRE_orbitphase_qflaginfo()
% calibrator_T: 1D array of internal cal target temperatures.
%     Should have the same shape as obs_ctime. [K]
% scipkt_ctime: 1D array of calibration-input science packet ctimes [s]
% ti: structure array containing various TIRS-related parameters
% TIRS_ID: integer, 1 = TIRS1, 2 = TIRS2
% has_bus_bitflags: (boolean) are bus status/event bitflags available in the
%     appropriate Cal_Artifacts arrays?
% Cal_Artifacts: structure array containing various calibration-related fields
%
% Outputs:
% obs_bitflags: 1D array of uint16 containing bitflags.
%     shape will match the input obs_ctime array.

% some threshold values - could eventually be moved into input
% config files, if needed.
post_PObS_timedelta = 360.0;  % [s]
post_pldoff_ewarm_timedelta = 10.0;  % [s]
post_pldoff_gtndT_timedelta = 25200.0;  % [s]
calseq_timedelta_thresh1 = 240.0;  % [s]
calseq_timedelta_thresh2 = 480.0;  % [s]
ROICintg_timedelta = -0.5*ti.ROIC_tau;  % [s]

if TIRS_ID == 1  % TIRS1
   cal_T_diff_threshold = 1.5;  % [K]
else  % TIRS2
   cal_T_diff_threshold = 1.5;  % [K]
end

% Initialize bitflag array:
obs_bitflags = zeros(length(obs_ctime), 1, 'uint16');

% Only if bus status and event bitflags are available:
if has_bus_bitflags
   % A mask for scipkt frames within this granule's time range, and a mask for
   %  obstgt frames within this granule's time range:
   granule_mask = ...
           (scipkt_ctime+ROICintg_timedelta > granule_beg_ts-0.001) & ...
           (scipkt_ctime+ROICintg_timedelta < granule_end_ts+0.001);
   granule_obstgt_mask = granule_mask & ...
                          (Cal_Artifacts.scipkt_frame_type == ti.OBSTGT_NOM);

   % Extract info from bus status bitflags (obstgt frames within granule):

   tmp = Cal_Artifacts.scipkt_bus_status_bitflags(granule_obstgt_mask);

   ATT_DET_invalid = logical(bitget(tmp, 1+1, 'uint8'));
   bus_tlm_gap = logical(bitget(tmp, 3+1, 'uint8'));

   % Extract info from bus status bitflags (within and before granule,
   %  any frame category):

     % Determine payload power-off period end ctimes:
   payload_powered_off = cast( ...
            bitget(Cal_Artifacts.scipkt_bus_status_bitflags, 4+1, 'uint8'), ...
                              'int32');
   end_of_pld_power_off_ctime = scipkt_ctime( ...
                                        find(diff(payload_powered_off) < 0)+1);

   % Extract info from bus events bitflags (obstgt frames within granule):

   tmp = Cal_Artifacts.scipkt_bus_events_bitflags(granule_obstgt_mask);

   unknown_bus_slew = logical(bitget(tmp, 6+1, 'uint16'));
   MSA_bus_slew = logical(bitget(tmp, 11+1, 'uint16'));

   % Extract info from bus events bitflags (within and before granule,
   %  any frame category):

     % Determine PObS (payload-on-but-safed; e.g., during/near downlink
     %  contacts and GPS zenith pointing) period end ctimes:
   payload_on_but_safed = cast( ...
            bitget(Cal_Artifacts.scipkt_bus_events_bitflags, 5+1, 'uint16'), ...
                               'int32');
   end_of_PObS_ctime = scipkt_ctime(find(diff(payload_on_but_safed) < 0)+1);

     % Determine eclipse entrance and exit ctimes:
   in_eclipse = cast( ...
            bitget(Cal_Artifacts.scipkt_bus_events_bitflags, 2+1, 'uint16'), ...
                               'int32');
   eclipse_entrance_ctime = scipkt_ctime(find(diff(in_eclipse) > 0)+1);
   eclipse_exit_ctime = scipkt_ctime(find(diff(in_eclipse) < 0));

     % (b0) Thermal transient after PObS (payload-on-but-safed) period (e.g.,
     %  during/near downlink contacts and GPS zenith pointing)
   post_PObS_flag_val = uint16(2.^0);  % (b0)
   for ipp=1:length(end_of_PObS_ctime)
      msk = (obs_ctime >= end_of_PObS_ctime(ipp)) & ...
            (obs_ctime < end_of_PObS_ctime(ipp)+post_PObS_timedelta);
      if any(msk)
         obs_bitflags(msk) = bitor(obs_bitflags(msk), post_PObS_flag_val);
      end
   end

   % (b6) bus telemetry indicates that spacecraft attitude determination is
   %  invalid
   hbt_attdet_invalid_flag_val = uint16(2.^6);  % (b6)
   msk = ATT_DET_invalid & ~bus_tlm_gap;
   if any(msk)
      obs_bitflags(msk) = bitor(obs_bitflags(msk), hbt_attdet_invalid_flag_val);
   end

   % (b7) lack of spacecraft attitude information due to a bus telemetry gap
   no_ia_bustlm_flag_val = uint16(2.^7);  % (b7)
   msk = bus_tlm_gap;
   if any(msk)
      obs_bitflags(msk) = bitor(obs_bitflags(msk), no_ia_bustlm_flag_val);
   end

   % (b8) during an unknown type of bus slew
   unknown_bus_slew_flag_val = uint16(2.^8);  % (b8)
   msk = unknown_bus_slew;
   if any(msk)
      obs_bitflags(msk) = bitor(obs_bitflags(msk), unknown_bus_slew_flag_val);
   end

   % (b9) during a modeled-sun-avoidance type of bus slew
   MSA_bus_slew_flag_val = uint16(2.^9);  % (b9)
   msk = MSA_bus_slew;
   if any(msk)
      obs_bitflags(msk) = bitor(obs_bitflags(msk), MSA_bus_slew_flag_val);
   end

   % (b10) Electronics warm-up period after powering instrument on
   % (b3) Greater-than-normal temperature change within this orbit (during the
   %  thermal transient after powering instrument back on)
   instr_ewarm_flag_val = uint16(2.^10);  % (b10)
   instr_temp_flag_val = uint16(2.^3);  % (b3)
   for ipp=1:length(end_of_pld_power_off_ctime)
      msk = (obs_ctime >= end_of_pld_power_off_ctime(ipp)) & ...
            (obs_ctime < end_of_pld_power_off_ctime(ipp)+ ...
             post_pldoff_ewarm_timedelta);
      if any(msk)
         obs_bitflags(msk) = bitor(obs_bitflags(msk), instr_ewarm_flag_val);
      end

      msk = (obs_ctime >= end_of_pld_power_off_ctime(ipp)) & ...
            (obs_ctime < end_of_pld_power_off_ctime(ipp)+ ...
             post_pldoff_gtndT_timedelta);
      if any(msk)
         min_cal_T = min(calibrator_T, [], 'all');  % Calculated over the
         max_cal_T = max(calibrator_T, [], 'all');  %  entire granule
           % Check if the calibrator T extrema were calculated over at least ~1
           %  orbit; if not, flag the obs even if cal_T_diff suggests not to:
         t_check = granule_end_ts-end_of_pld_power_off_ctime(ipp) < 5600.;
         if (max_cal_T-min_cal_T >= cal_T_diff_threshold) | t_check
            obs_bitflags(msk) = bitor(obs_bitflags(msk), instr_temp_flag_val);
         end
      end
   end

   % (b1) Small thermal/radiometric perturbation present (e.g., near eclipse
   %       exit)
   % (b2) Large thermal/radiometric perturbation present (e.g., near eclipse
   %       entrance)
   orbphase_Ee_flag_val = [uint16(2.^1), uint16(2.^2)];  % (b1), (b2)
   est_orbit_period = granule_end_ts-granule_beg_ts;  % [s]

   ophref_ct = [granule_beg_ts-est_orbit_period, granule_beg_ts, ...
                granule_end_ts, granule_end_ts+est_orbit_period];
   ophref = [0., 360., 720., 1080.];
   orbphase_tmp = orbphase+360.;

   Een_ct = eclipse_entrance_ctime(find(eclipse_entrance_ctime >= ...
                              granule_beg_ts-est_orbit_period*0.4, 3, 'first'));
   Een_oph = interp1(ophref_ct, ophref, Een_ct, 'linear', 'extrap');
   Eex_ct = eclipse_exit_ctime(find(eclipse_exit_ctime < ...
                              granule_end_ts+est_orbit_period*0.4, 3, 'last'));
   Eex_oph = interp1(ophref_ct, ophref, Eex_ct, 'linear', 'extrap');

   for n=1:length(orbphase_qinfo.level)
      if (orbphase_qinfo.level(n) == 0)
         continue
      end

      if orbphase_qinfo.ref(n) == 11
         for i=1:length(Een_oph)
            msk = (orbphase_tmp > ...
                   Een_oph(i)+orbphase_qinfo.orbitphase_offsets(n,1)) & ...
                  (orbphase_tmp < ...
                   Een_oph(i)+orbphase_qinfo.orbitphase_offsets(n,2));
            obs_bitflags(msk) = bitor(obs_bitflags(msk), ...
                                 orbphase_Ee_flag_val(orbphase_qinfo.level(n)));
         end
      elseif orbphase_qinfo.ref(n) == 12
         for i=1:length(Eex_oph)
            msk = (orbphase_tmp > ...
                   Eex_oph(i)+orbphase_qinfo.orbitphase_offsets(n,1)) & ...
                  (orbphase_tmp < ...
                   Eex_oph(i)+orbphase_qinfo.orbitphase_offsets(n,2));
            obs_bitflags(msk) = bitor(obs_bitflags(msk), ...
                                 orbphase_Ee_flag_val(orbphase_qinfo.level(n)));
         end
      end
   end
end

% (b4) Moderate time interval between observation and the nearest calibration
%       sequence (due to corrupted or missing calibration sequence(s))
% (b5) Large time interval between observation and the nearest calibration
%       sequence (due to corrupted or missing calibration sequence(s))
calseq_time_flag_val1 = uint16(2.^4);  % (b4)
calseq_time_flag_val2 = uint16(2.^5);  % (b5)
time_since_calseq = zeros(size(obs_ctime));
for t=1:length(obs_ctime)
   time_since_calseq(t) = min(abs(obs_ctime(t)-calseq_ctime), [], 'all');
end
msk1 = time_since_calseq > calseq_timedelta_thresh1;
if any(msk1)
   obs_bitflags(msk1) = bitor(obs_bitflags(msk1), calseq_time_flag_val1);
   msk2 = time_since_calseq > calseq_timedelta_thresh2;
   if any(msk2)
      obs_bitflags(msk2) = bitor(obs_bitflags(msk2), calseq_time_flag_val2);
   end
end

end
