"""
From one or more Level-0 payload telemetry file(s), categorize all science
data frames within, determining frame types and the locations of usable
calibration sequences.  Write the resulting information to a NetCDF-format file.

This program requires Python version 3.6 or later, and is importable as a
Python module.
"""

  # From the Python standard library:
import os
import datetime
from functools import reduce
from itertools import chain

  # From other external Python packages:
import numpy as np
import netCDF4
import scipy.ndimage

  # Custom utilities:
import PREFIRE_tools.TIRS.TIRS_packet as TIRSp
from PREFIRE_tools.utils.filesys import mkdir_p
from PREFIRE_tools.utils.time import ctime_to_UTC_DT
from PREFIRE_tools.utils.bitflags import apply_bit_to_bitflags_v, get_bit_from_bitflags_v
from PREFIRE_tools.utils.unit_conv import hour_to_s
from PREFIRE_PRD_GEN.file_creation import write_data_fromspec
import PREFIRE_L1.filepaths as L1_fpaths


debug = False

bus_vars_to_read = ["ctime", "REFS_is_valid", "REFS_ESM_is_valid",
                    "ATT_DET_is_valid",
                    "SC_in_moon_umbra", "q_scbody_wrt_ECI_1",
                    "q_scbody_wrt_ECI_2", "q_scbody_wrt_ECI_3",
                    "q_scbody_wrt_ECI_4",
                    "pri_cmd_ref_dir", "pri_cmd_vec_wrt_SBF_1",
                    "cmd_q_scbody_wrt_ECI_1", "cmd_q_scbody_wrt_ECI_2",
                    "cmd_q_scbody_wrt_ECI_3", "cmd_q_scbody_wrt_ECI_4",
                    "radio_is_tx", "radio_rx_lock", "MSA_action_occurring",
                    "payload_is_on"]

orbit_vars_to_read = ["ctime", "SC_in_Earth_umbra", "velocity_wrt_ECI_3"]

# bus slew to contact attitude (1), bus slew to LVLH (2),
# bus slew to safe attitude (3), bus slew to safe attitude (3),
# bus slew to contact attitude (1), bus slew to LVLH (2)
# bus slew to other/unknown attitude (0), bus slew to LVLH (2)
pri_cmd_chg_t_ref_d = {(13, 9): 1, (9, 13): 2,
                       (13, 1): 3, (9, 1): 3,
                       (1, 9): 1, (1, 13): 2,
                       (13, 3): 0, (3, 13): 2}


#--------------------------------------------------------------------------
def general_tlm_file_info(input_fpaths, tlm_type, start_ctime=None,
                          end_ctime=None):
    """Extract/determine general info about a set of prepped telemetry files."""

    # Expects input NetCDF-format prepped telemetry files in ascending
    #  chronological order (both lexicographically by filename, and internally).
    #
    # start_ctime : (optional) float64 or None
    #      ctime at which to start reading data.  All input frames that have a
    #       ctime less than this will be ignored.
    # end_ctime : (optional) float64 or None
    #      ctime at which to stop reading data.  All input frames that have a
    #       ctime greater than this will be ignored.

    tmp_sum = sum([start_ctime is not None, end_ctime is not None])
    if tmp_sum > 0.99 and tmp_sum < 1.01:
        raise ValueError("start_ctime and end_ctime must either both be None, "
                         "or both be floating point values.")

    subset_in_ctime_status = [False, False, False]
    subset_in_ctime_status[0] = ( tmp_sum > 1.5 )

    if tlm_type == "payload":
        atrack_dim_mnk = "atrack_scipkt"
        dgroup_mnk = "TIRS_L0_Data"
    elif tlm_type == "bus":
        atrack_dim_mnk = "atrack"
        dgroup_mnk = "SCbus_L0_Data"
    else:
        atrack_dim_mnk = "atrack"
        dgroup_mnk = "Orbit_L0_Data"

    ds_l = []
    fn_l = []
    SC_ID_set = set()
    sensor_ID_set = set()
    n_atrack_tl = []
    for input_fpath in input_fpaths:
        ds = netCDF4.Dataset(input_fpath)

        try:
            n_atrack = ds.dimensions[atrack_dim_mnk].size
        except:
            n_atrack = ds.groups[dgroup_mnk].dimensions[atrack_dim_mnk].size
        ict_cov_l = [0, n_atrack]  # Default, NumPy indexing convention

        if subset_in_ctime_status[0]:
            ct_cov_t = (ds.ctime_coverage_start_s, ds.ctime_coverage_end_s)
            if (ct_cov_t[0] > end_ctime or ct_cov_t[1] < start_ctime):
                continue  # This input file is not needed

            if (start_ctime >= ct_cov_t[0] and start_ctime <= ct_cov_t[1]):
                ctime = ds.groups[dgroup_mnk].variables["ctime"]
                ict_cov_l[0] = np.searchsorted(ctime, start_ctime, side="left")
                subset_in_ctime_status[1] = True
            if (end_ctime >= ct_cov_t[0] and end_ctime <= ct_cov_t[1]):
                try:
                    t = len(ctime)  # Just to check if ctime is instantiated
                except:
                    ctime = ds.groups[dgroup_mnk].variables["ctime"]
                ict_cov_l[1] = np.searchsorted(ctime, end_ctime, side="right")
                subset_in_ctime_status[2] = True

        ds_l.append(ds)
        fn_l.append(os.path.basename(input_fpath))
        n_atrack_tl.append((ict_cov_l[1]-ict_cov_l[0], ict_cov_l[0],
                            ict_cov_l[1]))

        orbit_sim_version = ds.orbit_sim_version
        SRF_NEdR_version = ds.SRF_NEdR_version
        SC_ID_set.add(ds.spacecraft_ID)
        sensor_ID_set.add(ds.sensor_ID)

    if subset_in_ctime_status[0]:
        chk = sum(subset_in_ctime_status[1:3])
        if chk < 2:
            ct_s = sorted([ds.ctime_coverage_start_s for ds in ds_l])
            ct_e = sorted([ds.ctime_coverage_end_s for ds in ds_l])
            if chk == 0:
                rsn = "both the start and end of the range"
            elif not subset_in_ctime_status[1]:
                rsn = "the start of the range"
            elif not subset_in_ctime_status[2]:
                rsn = "the end of the range"
            rsn += (" (time range of supplied input data: ctime {} to {})"
                                                    .format(ct_s[0], ct_e[-1]))
            msg = ("Unable to extract a time subset (ctime {} to {}) of data "
                   "group {}, likely due to insufficient input data supplied "
                   "at {}.".format(start_ctime, end_ctime, dgroup_mnk, rsn))
            if ct_s[0]-start_ctime > 2. or end_ctime-ct_e[-1] > 2.:
                # Outside a nominal 2 s margin of error -- raise an exception:
                raise RuntimeError(msg)

    if len(sensor_ID_set) == 1:
        sensor_ID = min(sensor_ID_set)
        if len(sensor_ID.strip()) > 0:
            i_tmp = int(sensor_ID[5])
            sensor_IDval, sensor_IDidx = (i_tmp, i_tmp-1)
        else:
            sensor_IDval, sensor_IDidx = (None, None)
    else:
        raise ValueError("'sensor_ID' is not identical for all input "
                         "files ({})".format(','.join(input_fpaths)))

    if len(SC_ID_set) == 1:
        SC_ID = min(SC_ID_set)
    else:
        raise ValueError("'spacecraft_ID' is not identical for all input "
                         "files ({})".format(','.join(input_fpaths)))

    n_atrack_total = sum([x[0] for x in n_atrack_tl])

    if tlm_type == "payload":
        try:
            n_xtrack = ds_l[0].dimensions["xtrack"].size
        except:
            try:
                n_xtrack = ds_l[0].groups[dgroup_mnk].dimensions["xtrack"].size
            except:
                n_xtrack = None  # This is prototlm
    else:
        n_xtrack = None

    return {"ds_l": ds_l, "fn_l": fn_l, "n_atrack_total": n_atrack_total,
            "n_atrack_tl": n_atrack_tl,
            "orbit_sim_version": orbit_sim_version,
            "SRF_NEdR_version": SRF_NEdR_version, "SC_ID": SC_ID,
            "sensor_ID": sensor_ID, "sensor_IDval": sensor_IDval,
            "sensor_IDidx": sensor_IDidx, "n_xtrack": n_xtrack}


#--------------------------------------------------------------------------
def _find_via_encoder(encoder, frametype_encoder_vals, use_xor_neq=False,
                      output_slices=False):
    if use_xor_neq:
        msk = reduce(np.logical_xor,
                     [encoder != x for x in frametype_encoder_vals])
    else:
        msk = reduce(np.logical_or,
                     [encoder == x for x in frametype_encoder_vals])
    some_possible_fr = np.any(msk)  # Do any frames match?
    if some_possible_fr:
        lbls, _ = scipy.ndimage.label(msk)
        slices = [s[0] for s in scipy.ndimage.find_objects(lbls)]
        fr_ind = np.r_[tuple(slices)]
    else:
        fr_ind = None
        slices = None

    if output_slices:
        return (some_possible_fr, slices)
    else:
        return (some_possible_fr, fr_ind)


#--------------------------------------------------------------------------
def _find_and_vet_calseqpart_edges_A(fr_ctime, min_n_frames_per_block):
    """ """

    # A space/caltgt-stare edge is indicated by > 25 s between
    #  space/caltgt frames, where that threshold more than accomodates
    #    * motor rotation time on both sides of the space/caltgt stare
    #    * up to 24 space/caltgt frames (up to ~16.82 s) within the stare
    # This version DOES require immediately-adjacent space and
    #  caltgt stares.

    ia_interior_block_edges = np.nonzero(np.diff(fr_ctime) > 25.)[0]
    if len(ia_interior_block_edges) > 0:
        ia_block_ld_edges0 = [0]  # leading edge, by default
        ia_block_tr_edges0 = [ia_interior_block_edges[0]]
        ia_block_ld_edges0.extend(ia_interior_block_edges[:-1]+1)  # Inner
        ia_block_tr_edges0.extend(ia_interior_block_edges[1:])  # Inner
        ia_block_ld_edges0.append(ia_interior_block_edges[-1]+1)
        ia_block_tr_edges0.append(len(fr_ctime)-1)  # trailing edge
    else:  # No interior block edges
        ia_block_ld_edges0 = [0]  # leading edge
        ia_block_tr_edges0 = [len(fr_ctime)-1]  # trailing edge
    ia_block_ld_edges = np.array(ia_block_ld_edges0)
    ia_block_tr_edges = np.array(ia_block_tr_edges0)

    # Reject blocks containing < 'min_n_frames_per_block' possible
    #  space/caltgt frames:
    ik_good_block_edges = np.nonzero(ia_block_tr_edges-ia_block_ld_edges+1 >=
                                                      min_n_frames_per_block)[0]
    if len(ik_good_block_edges) > 0:  # Some good blocks were found
        ic_spc_ld_edges = ia_block_ld_edges[ik_good_block_edges]
        ic_spc_tr_edges = ia_block_tr_edges[ik_good_block_edges]
        return (ic_spc_ld_edges, ic_spc_tr_edges)
    else:
        return None


def _find_and_vet_calseqpart_edges_B(fr_ctime, min_n_frames_per_block):
    """ """

    # A space/caltgt-stare edge is indicated by > 25 s between
    #  space/caltgt frames, where that threshold more than accomodates
    #    * motor rotation time on both sides of the space/caltgt stare
    #    * up to 24 space/caltgt frames (up to ~16.82 s) within the stare
    # This version DOES NOT require immediately-adjacent space and
    #  caltgt stares."""

    ia_interior_block_edges = np.nonzero(np.diff(fr_ctime) > 25.)[0]
    if len(ia_interior_block_edges) > 0:
        ia_block_ld_edges0 = [0]  # leading edge, by default
        ia_block_tr_edges0 = [ia_interior_block_edges[0]]
        ia_block_ld_edges0.extend(ia_interior_block_edges[:-1]+1)  # Inner
        ia_block_tr_edges0.extend(ia_interior_block_edges[1:])  # Inner
        ia_block_ld_edges0.append(ia_interior_block_edges[-1]+1)
        ia_block_tr_edges0.append(len(fr_ctime)-1)  # trailing edge
    else:  # No interior block edges
        ia_block_ld_edges0 = [0]  # leading edge
        ia_block_tr_edges0 = [len(fr_ctime)-1]  # trailing edge
    ia_block_ld_edges = np.array(ia_block_ld_edges0)
    ia_block_tr_edges = np.array(ia_block_tr_edges0)

    # Reject blocks containing < 'min_n_frames_per_block' possible
    #  space/caltgt frames:
    ik_good_block_edges = np.nonzero(ia_block_tr_edges-ia_block_ld_edges+1 >=
                                                      min_n_frames_per_block)[0]
    if len(ik_good_block_edges) > 0:  # Some good blocks were found
        ic_spc_ld_edges = ia_block_ld_edges[ik_good_block_edges]
        ic_spc_tr_edges = ia_block_tr_edges[ik_good_block_edges]
        return (ic_spc_ld_edges, ic_spc_tr_edges)
    else:
        return None


#--------------------------------------------------------------------------
def _find_invalid_calseqpart_fr_A(tgt_ltID, frametype, ind0, ctime, DN_c12_s4,
                                  DN_c12_s6):
    """Find invalid frames in a portion (space or caltgt) of the calibration
       sequence.  This version DOES require immediately-adjacent space and
       caltgt stares."""
    ctime0 = ctime[ind0]
    w_frametype = frametype.copy()

    # Find space/caltgt-stare edges, rejecting any blocks containing < 5
    #  possible space/caltgt frames:
    edges_t = _find_and_vet_calseqpart_edges_A(ctime0, 5)
    if edges_t is None:
        return (w_frametype, edges_t, None)  # No good blocks; done here

    ic_ld_edges, ic_tr_edges = edges_t
    
      # Deem any frames culled at this step to be skew, and record them as
      #  such in the 'w_frametype' array:
    w_frametype[ind0] = TIRSp.SKEW_VIA_BLKCNT1
    ind_msk = np.r_[tuple([slice(i_ld, i_tr+1) for i_ld, i_tr in
                                               zip(ic_ld_edges, ic_tr_edges)])]
    w_frametype[ind0[0][ind_msk]] = tgt_ltID

    DN_c12_s4_0 = DN_c12_s4[ind0]
    DN_c12_s6_0 = DN_c12_s6[ind0]

    block_slice_l = [slice(ic_ld_edges[i], ic_tr_edges[i]+1)
                                              for i in range(len(ic_ld_edges))]
      # Compute (DN - median_DN_of_block) for every block of frames:
    diff_DN_c12_s4 = np.array(list(chain.from_iterable(
                [np.abs(DN_c12_s4_0[s]-np.median(DN_c12_s4_0[s])) for s in
                                                              block_slice_l])))
    diff_DN_c12_s6 = np.array(list(chain.from_iterable(
                [np.abs(DN_c12_s6_0[s]-np.median(DN_c12_s6_0[s])) for s in
                                                              block_slice_l])))

      # Frames with values of (DN - median_DN_of_block) > 9 are deemed to
      #  be skew -- record them as such in the 'w_frametype' array:
    block_ind = np.array(np.r_[tuple(block_slice_l)], dtype="int32")
    ind_msk = np.nonzero((diff_DN_c12_s4 > 9) & (diff_DN_c12_s6 > 9))
    w_frametype[ind0[0][block_ind[ind_msk]]] = TIRSp.SKEW_VIA_DN

    # Re-inventory space/caltgt-stare edges, rejecting any blocks
    #  containing < 3 possible space/caltgt frames:
    ind1 = np.nonzero(w_frametype == tgt_ltID)
    edges_t = _find_and_vet_calseqpart_edges_A(ctime[ind1], 3)
    if edges_t is None:
        return (w_frametype, edges_t, None)  # No good blocks; done here

      # Deem any frames culled at this step to be skew, and record them as
      #  such in the 'w_frametype' array:
    w_frametype[ind1] = TIRSp.SKEW_VIA_BLKCNT2
    ind_msk = np.r_[tuple([slice(i_ld, i_tr+1) for i_ld, i_tr in
                                                 zip(edges_t[0], edges_t[1])])]
    w_frametype[ind1[0][ind_msk]] = tgt_ltID

      # Construct edge and block-frame indices relative to the 'w_frametype'
      #  array:
    abs_edges_t = (np.array([ind1[0][i] for i in edges_t[0]]),
                   np.array([ind1[0][i] for i in edges_t[1]]))
    abs_block_ind = [list(range(i_ld, i_tr+1)) for i_ld, i_tr in
                                           zip(abs_edges_t[0], abs_edges_t[1])]

    # A frame within 'ROIC_tau' seconds directly after the last caltgt
    #  frame (or any skew frame adjacent to that) should probably be
    #  considered skew, out of an abundance of caution:
    #   (note that due to the caltgt frame block size constraints, issues
    #    regarding being too close to the array beginning do not need to be
    #    treated here)
    if tgt_ltID == TIRSp.CALTGT_POSS:
        n_frames_all = len(ctime)
        for i_caltgtb in range(len(abs_edges_t[1])):
            it = abs_block_ind[i_caltgtb][-1]
            if it < n_frames_all-1:  # if not too close to the array end
                it += 1
                if w_frametype[it]//10 == 3:  # Skew frame
                    it += 1
                if it <= n_frames_all-1:  # if index still within the array
                    if (ctime[it]-ctime[it-1]) < TIRSp.ROIC_tau+0.1:
                        w_frametype[it] = TIRSp.SKEW_VIA_DISP

    return (w_frametype, abs_edges_t, abs_block_ind)


#--------------------------------------------------------------------------
def _find_invalid_calseqpart_fr_B(tgt_ltID, frametype, ind0, ctime, DN_c12_s4,
                                  DN_c12_s6, use_min_DN=False):
    """Find invalid frames in a portion (space or caltgt) of the calibration
       sequence. This version DOES NOT require immediately-adjacent space and
       caltgt stares."""
    ctime0 = ctime[ind0]
    w_frametype = frametype.copy()

    # Find space/caltgt-stare edges, rejecting any blocks containing < 4
    #  possible space/caltgt frames:
    edges_t = _find_and_vet_calseqpart_edges_B(ctime0, 4)
    if edges_t is None:
        return (w_frametype, edges_t, None)  # No good blocks; done here

    ic_ld_edges, ic_tr_edges = edges_t
    
      # Deem any frames culled at this step to be skew, and record them as
      #  such in the 'w_frametype' array:
    w_frametype[ind0] = TIRSp.SKEW_VIA_BLKCNT1
    ind_msk = np.r_[tuple([slice(i_ld, i_tr+1) for i_ld, i_tr in
                                               zip(ic_ld_edges, ic_tr_edges)])]
    w_frametype[ind0[0][ind_msk]] = tgt_ltID

    DN_c12_s4_0 = DN_c12_s4[ind0]
    DN_c12_s6_0 = DN_c12_s6[ind0]

    block_slice_l = [slice(ic_ld_edges[i], ic_tr_edges[i]+1)
                                              for i in range(len(ic_ld_edges))]
    if use_min_DN:
        # Compute (DN - min_DN_of_block) for every block of frames:
        diff_DN_c12_s4 = np.array(list(chain.from_iterable(
                       [np.abs(DN_c12_s4_0[s]-np.amin(DN_c12_s4_0[s])) for s
                                                           in block_slice_l])))
        diff_DN_c12_s6 = np.array(list(chain.from_iterable(
                       [np.abs(DN_c12_s6_0[s]-np.amin(DN_c12_s6_0[s])) for s
                                                           in block_slice_l])))
    else:
        # Compute (DN - median_DN_of_block) for every block of frames:
        diff_DN_c12_s4 = np.array(list(chain.from_iterable(
                       [np.abs(DN_c12_s4_0[s]-np.median(DN_c12_s4_0[s])) for s
                                                           in block_slice_l])))
        diff_DN_c12_s6 = np.array(list(chain.from_iterable(
                       [np.abs(DN_c12_s6_0[s]-np.median(DN_c12_s6_0[s])) for s
                                                           in block_slice_l])))

      # Frames with values of (DN - median_DN_of_block) > 11 are deemed to
      #  be skew -- record them as such in the 'w_frametype' array:
    dDN_th = 11
    block_ind = np.array(np.r_[tuple(block_slice_l)], dtype="int32")
    ind_msk = np.nonzero((diff_DN_c12_s4 > dDN_th) & (diff_DN_c12_s6 > dDN_th))
    w_frametype[ind0[0][block_ind[ind_msk]]] = TIRSp.SKEW_VIA_DN

    # Re-inventory space/caltgt-stare edges, rejecting any blocks
    #  containing < 3 possible space/caltgt frames:
    ind1 = np.nonzero(w_frametype == tgt_ltID)
    edges_t = _find_and_vet_calseqpart_edges_A(ctime[ind1], 3)
    if edges_t is None:
        return (w_frametype, edges_t, None)  # No good blocks; done here

      # Deem any frames culled at this step to be skew, and record them as
      #  such in the 'w_frametype' array:
    w_frametype[ind1] = TIRSp.SKEW_VIA_BLKCNT2
    ind_msk = np.r_[tuple([slice(i_ld, i_tr+1) for i_ld, i_tr in
                                                 zip(edges_t[0], edges_t[1])])]
    w_frametype[ind1[0][ind_msk]] = tgt_ltID

      # Construct edge and block-frame indices relative to the 'w_frametype'
      #  array:
    abs_edges_t = (np.array([ind1[0][i] for i in edges_t[0]]),
                   np.array([ind1[0][i] for i in edges_t[1]]))
    abs_block_ind = [list(range(i_ld, i_tr+1)) for i_ld, i_tr in
                                           zip(abs_edges_t[0], abs_edges_t[1])]

    # Record "drive" frames and unused calseq frames as being SKEW_VIA_DISP:
    if tgt_ltID == TIRSp.CALTGT_POSS:
        n_frames_all = len(ctime)
          # Very conservative (excessive in most situations?) duration estimates
          #  for the first (space stare triplet + drives) and
          #  second (caltgt stare and drive) parts of the calseq:
        t_1pcs = 18.7-3.7  # [s]
        t_2pcs = 5.2  # [s]

        for i_caltgtb in range(len(abs_edges_t[1])):
            itb = abs_block_ind[i_caltgtb][0]
               # Search for a "drive" frame between the space-stare triplet and
               #  the first caltgt stare frame -- if one is not identifiable
               #  within 4 frames, use the frame before the first caltgt stare
               #  frame):
            tmp_th = DN_c12_s4[itb]-(dDN_th+1)
            tmp = [i for i in range(max(0, itb-4), itb) if
                      ((DN_c12_s4[i] <= tmp_th) and (w_frametype[i]//10 != 2))]
            if len(tmp) > 0:
                isc = tmp[-1]
            else:
                isc = max(0, itb-1)
            ctime_isc = ctime[isc]
               # Step backward from the identified frame to the potentially
               #  excessive estimate of the drive start to the first space stare
               #  of this calseq, marking all OBSTGT_NOM frames encountered as
               #  SKEW_VIA_DISP frames:
            start_of_cseq_ct = ctime_isc-t_1pcs
            it = isc
            while (ctime[it] >= start_of_cseq_ct):
                if w_frametype[it] == TIRSp.OBSTGT_NOM:
                    w_frametype[it] = TIRSp.SKEW_VIA_DISP
                it -= 1
                if it < 0:
                    break

            # Step forward from the identified frame to the potentially
            #  excessive estimate of the end of the caltgt stare and drive part
            #  of the calseq, marking all OBSTGT_NOM frames encountered as
            #  SKEW_VIA_DISP frames:
            end_of_cseq_ct = ctime_isc+t_2pcs
            ite = abs_block_ind[i_caltgtb][-1]
            it = isc+1
            while (ctime[it] <= end_of_cseq_ct):
                if w_frametype[it] == TIRSp.OBSTGT_NOM:
                    w_frametype[it] = TIRSp.SKEW_VIA_DISP
                it += 1
                if it > n_frames_all-1:
                    break

    return (w_frametype, abs_edges_t, abs_block_ind)


#--------------------------------------------------------------------------
def sift_SAT2_data(frametype, space_poss_slices, bus_status_bfs, bus_events_bfs,
                   HS_slices, ctime, encoder_pos, DN_c12_s4, DN_c12_s6, TIRS_T):
    """Sift through the available SAT2 data, identifying/refining calibration
       sequence frames, etc."""

    n_calseq_stare_fr = [4, 9]
    calseq_stare_max_t_duration = n_calseq_stare_fr[1]*TIRSp.ROIC_tau*1.02
    space_sample_quantile = 0.50
    deltaDN_calseq_th = 1300.

    mirror_aostart_ctime = []
    mirror_ang_offset = []

    for s in space_poss_slices:
        s0 = s[0]
        n_fr = s0.stop-s0.start

        if (n_fr >= n_calseq_stare_fr[0] and n_fr <= n_calseq_stare_fr[1] and
            ctime[s0.stop-1]-ctime[s0.start] < calseq_stare_max_t_duration):
            # This cal-sequence space stare candidate is of an acceptable frame
            #  length and time duration; vet it further:
            passed, frametype = vet_SAT2_sps(s0, s, DN_c12_s4, DN_c12_s6,
                         deltaDN_calseq_th, space_sample_quantile, frametype,
                         ctime, calseq_stare_max_t_duration, n_calseq_stare_fr)

    # Determine invalid or unusable apparent/possible space frames:
    ind_space0 = np.nonzero(np.abs(frametype) == TIRSp.SPACE_POSS)
    if len(ind_space0[0]) > 0:
        frametype, space_edge_t, space_block_ind = (
                    _find_invalid_calseqpart_fr_A(TIRSp.SPACE_POSS, frametype,
                                                  ind_space0, ctime, DN_c12_s4,
                                                  DN_c12_s6))
        frametype[frametype == -TIRSp.SPACE_POSS] = TIRSp.SPACE_POSS
    else:
        space_edge_t, space_block_ind = (None, None)

    # Determine invalid or unusable apparent/possible caltgt frames:
    ind_caltgt0 = np.nonzero(np.abs(frametype) == TIRSp.CALTGT_POSS)
    if len(ind_caltgt0[0]) > 0:
        frametype, caltgt_edge_t, caltgt_block_ind = (
                   _find_invalid_calseqpart_fr_A(TIRSp.CALTGT_POSS, frametype,
                                                 ind_caltgt0, ctime, DN_c12_s4,
                                                 DN_c12_s6))
        frametype[frametype == -TIRSp.CALTGT_POSS] = TIRSp.CALTGT_POSS
    else:
        caltgt_edge_t, caltgt_block_ind = (None, None)

    if len(mirror_aostart_ctime) == 0:
        mirror_aostart_ctime, mirror_ang_offset = (None, None)
        
    return (frametype, (space_edge_t, space_block_ind),
            (caltgt_edge_t, caltgt_block_ind),
            (mirror_aostart_ctime, mirror_ang_offset))


#--------------------------------------------------------------------------
def vet_SAT2_sps(s0, s, DN_c12_s4, DN_c12_s6, deltaDN_calseq_th,
                 space_sample_quantile, frametype, ctime,
                 calseq_stare_max_t_duration, n_calseq_stare_fr):
    """Vet a possible SAT2 cal-sequence space stare candidate."""

    passed = False
    
    # Is there an immediately-enough-following caltgt stare candidate?
    idb = s0.stop+2
    ide = min(idb+3, len(frametype))  # Select <=3 frames; NumPy indexing
    if idb >= ide-1:
        return (passed, frametype)  # No suitable caltgt stare found
    if (any(frametype[idb:ide] != TIRSp.CALTGT_POSS) or
        ctime[ide-1]-ctime[idb] >= calseq_stare_max_t_duration):
        return (passed, frametype)  # No suitable caltgt stare found
    space_sample_DN_s4 = np.quantile(DN_c12_s4[s], space_sample_quantile)
    enough_dDN_s4 = ( np.median(DN_c12_s4[idb:ide])-space_sample_DN_s4 >
                                                            deltaDN_calseq_th )
    space_sample_DN_s6 = np.quantile(DN_c12_s6[s], space_sample_quantile)
    enough_dDN_s6 = ( np.median(DN_c12_s6[idb:ide])-space_sample_DN_s6 >
                                                            deltaDN_calseq_th )
    if enough_dDN_s4 or enough_dDN_s6:
        # This is a "good enough" cal-sequence space stare:
        passed = True
        frametype[s] = -TIRSp.SPACE_POSS

        # Determine any associated caltgt-stare frames:
        idb = s0.stop
        ide = min(idb+n_calseq_stare_fr[1], len(frametype))
        tmp_i = np.nonzero(frametype[idb:ide] == TIRSp.CALTGT_POSS)
        frametype[tmp_i[0][0]+idb:tmp_i[0][-1]+idb+1] = -TIRSp.CALTGT_POSS

    return (passed, frametype)


#--------------------------------------------------------------------------
def vet_SAT1_sps(s0, s, DN_c12_s4, DN_c12_s6, deltaDN_calseq_info,
                 space_sample_quantile, frametype, ctime,
                 calseq_sp_stare_max_t_duration, n_calseq_sp_stare_fr,
                 n_calseq_ca_stare_fr, TIRS_T, mirror_angle_offset):
    """Vet a possible SAT1 cal-sequence space stare candidate."""

    passed = False
    i_T_caltgt = slice(3, 6)
    i_Een_T_active = 1
    i_Een_T_passive = 0

    # Is there an immediately-enough-following caltgt stare candidate?

    space_sample_DN_s4 = np.quantile(DN_c12_s4[s], space_sample_quantile)
    space_sample_DN_s6 = np.quantile(DN_c12_s6[s], space_sample_quantile)

    n_s = s0.stop-s0.start  # Number of frames in space stare

       # Search for caltgt stare candidates in this index range:
    idb = s0.stop+1
    ide = min(idb+n_calseq_sp_stare_fr[1], len(ctime))  # NumPy indexing

    # AJM indicates that the linear fit (delta-DN with caltgt temperature)
    #  should hold "well enough" under normal orbital operating temperatures.
    # (However, all bets are off during periods when the caltgt T is not
    #  representative of the optical plane T.)
    eff_caltgt_T = np.mean(TIRS_T[idb:ide,i_T_caltgt], axis=1)

    Een_T_ratio = (TIRS_T[s0.start:ide,i_Een_T_active]/
                   TIRS_T[s0.start:ide,i_Een_T_passive])
    shift_factor = [0, 1, -1, 2, -2, 3, -3, 4, -4]

    if debug:
        if (ctime[s0.start] >= ctime[-1]-3.2*3600.-1500. and
            ctime[s0.start] <= ctime[-1]-1.6*3600.+1500.):
            print("ccc",s0.start,ctime[s0.start]-(ctime[-1]-3.2*3600.-1500.),
                  np.mean(Een_T_ratio))

    dDN_s4 = DN_c12_s4[idb:ide]-space_sample_DN_s4
    dDN_s4_th0 = (deltaDN_calseq_info["offset"][0]+
                  deltaDN_calseq_info["slope"][0]*eff_caltgt_T)
    dDN_s4_th_tol = dDN_s4_th0*0.02
    if np.mean(Een_T_ratio) > 1.00081:
        for i, shift in enumerate(shift_factor):
            dDN_s4_th = dDN_s4_th0+shift*dDN_s4_th_tol
            ii4 = np.nonzero((dDN_s4 >= dDN_s4_th-dDN_s4_th_tol) &
                             (dDN_s4 <= dDN_s4_th+dDN_s4_th_tol))
            if len(ii4[0]) >= 5:
                if i > 0:
                    print("nonzero sc4 DN threshold shift that appears to work",
                          i, shift, space_sample_DN_s4, len(ii4[0]))
                break
    else:
        ii4 = np.nonzero((dDN_s4 >= dDN_s4_th0-dDN_s4_th_tol) &
                         (dDN_s4 <= dDN_s4_th0+dDN_s4_th_tol))

    dDN_s6 = DN_c12_s6[idb:ide]-space_sample_DN_s6
    dDN_s6_th0 = (deltaDN_calseq_info["offset"][1]+
                 deltaDN_calseq_info["slope"][1]*eff_caltgt_T)
    dDN_s6_th_tol = dDN_s6_th0*0.02
    if np.mean(Een_T_ratio) > 1.00081:
        for i, shift in enumerate(shift_factor):
            dDN_s6_th = dDN_s6_th0+shift*dDN_s6_th_tol
            ii6 = np.nonzero((dDN_s6 >= dDN_s6_th-dDN_s6_th_tol) &
                             (dDN_s6 <= dDN_s6_th+dDN_s6_th_tol))
            if len(ii6[0]) >= 5:
                if i > 0:
                    print("nonzero sc6 DN threshold shift that appears to work",
                          i, shift, space_sample_DN_s6, len(ii6[0]))
                break
    else:
        ii6 = np.nonzero((dDN_s6 >= dDN_s6_th0-dDN_s6_th_tol) &
                         (dDN_s6 <= dDN_s6_th0+dDN_s6_th_tol))

    if len(ii4[0]) == 0 and len(ii6[0]) == 0:
        return (passed, frametype)

    # These seem to be "good enough" space stare and caltgt stare candidates:

    min_i, max_i = (len(ctime), -1)
    if len(ii4[0]) > 0:
        min_i, max_i = (min(min_i, ii4[0][0]), max(max_i, ii4[0][-1]))
    if len(ii6[0]) > 0:
        min_i, max_i = (min(min_i, ii6[0][0]), max(max_i, ii6[0][-1]))
    s_caltgt = slice(min_i+idb,
                     min(len(ctime), max_i+idb, min_i+idb+n_calseq_ca_stare_fr))
       # Check that the candidate caltgt frames are not too far apart in time:
    if (ctime[s_caltgt.stop-1]-ctime[s_caltgt.start] <
                                       TIRSp.ROIC_tau*(n_calseq_ca_stare_fr+2)):
        itb = s_caltgt.start
        # Search for a "drive" frame between the space-stare triplet and
        #  the first caltgt stare frame (within 4 frames of this):
        tmp = []
        ilow = max(0, itb-4)
        if ilow < itb:
            # Look for a high-amplitude positive DN difference beween frames.
            dDNtmp = np.diff(DN_c12_s4[ilow:itb+1])
            if mirror_angle_offset < 1.:
                DN_th_tmp = 100
            elif mirror_angle_offset < 1.1:
                DN_th_tmp = 50
            else:
                DN_th_tmp = 30
            tmp = [i+ilow for i in range(len(dDNtmp)) if dDNtmp[i] >= DN_th_tmp]
            if debug:
               print("dDN", mirror_angle_offset, np.min(dDNtmp), np.max(dDNtmp))
               if len(tmp) > 0:
                   print("ggg", np.min([DN_c12_s4[x+1]-DN_c12_s4[x]
                                        for x in tmp]))

        if len(tmp) > 0:
            # Identified a "drive" frame just before caltgt stare.
            passed = True
            frametype[s_caltgt] = TIRSp.CALTGT_POSS
            frametype[s] = TIRSp.SPACE_POSS

    return (passed, frametype)


#--------------------------------------------------------------------------
def sift_SAT1_data(frametype, space_poss_slices, bus_status_bfs, bus_events_bfs,
                   HS_slices, ctime, encoder_pos, DN_c12_s4, DN_c12_s6, TIRS_T):
    """Sift through the available SAT1 data, identifying/refining calibration
       sequence frames, etc."""

    # For calibration sequence space stare:
    n_calseq_sp_stare_fr = [4, 27]  # [min, max] frames
    calseq_sp_stare_max_t_duration = (
                             n_calseq_sp_stare_fr[1]*TIRSp.ROIC_tau*1.02)  # [s]
    space_sample_quantile = 0.20

    # For calibration sequence caltgt stare:
    n_calseq_ca_stare_fr = 6  # frames

    SASS_angspd = 0.5625  # [deg/s] (PHD 4000; 8x faster than A1 slow scan)
    SASS_ref_rel_ct = 18.94  # [s] from hard stop trailing edge (mean of SASS
                      #  space centers from 25 July to 24 August 2024; from NBM)

    # From just-before-launch analysis (and resulting data file
    #  TIRS1_prelaunch_RadCal_deltaDN.h5) by AJM for all channels and scenes:
    # "The deltaDN is the deltaDN (internal Cal - simulated space) during some
    #  number of calseq's during the Rad Cal test. The temperature of the
    #  internal cal is given by the calibT array. The calibT range for both
    #  TIRS tests was roughly 286 - 291 K. Within that range, the change in
    #  channel radiance is roughly linear - so I just did linear fits between
    #  the deltaDN and the calibT. The results of that are stored in the
    #  deltaDN_offset and deltaDN_slope arrays."
    #
    # The below offset & slope values are for channel 12, scenes 4 and 6:
    deltaDN_calseq_info = {"offset": (-3831.0137778099993, -3854.5890583911846),
                           "slope": (17.54756232496945, 17.54021982362612),
                           "T_bounds": (286.39447021, 290.92370605)}

    mirror_aostart_ctime = []
    mirror_ang_offset = []

    n_HS, n_sps = (len(HS_slices), len(space_poss_slices))
 
    # Case without any HS, or possible space stare slices before first HS:
    if n_HS > 0:
        i_HS = 0
        i_of_HS_edge = HS_slices[i_HS].stop-1
    else:
        i_of_HS_edge = 999999999
    i_sps = 0
    s0 = space_poss_slices[i_sps][0]
    while s0.start <= i_of_HS_edge:
        n_fr = s0.stop-s0.start
        if (n_fr >= n_calseq_sp_stare_fr[0] and
            n_fr <= n_calseq_sp_stare_fr[1] and
            ctime[s0.stop-1]-ctime[s0.start] < calseq_sp_stare_max_t_duration):
            # This cal-sequence space stare candidate is of an acceptable frame
            #  length and time duration; vet it further:
            passed, frametype = vet_SAT1_sps(s0, space_poss_slices[i_sps],
                                DN_c12_s4, DN_c12_s6, deltaDN_calseq_info,
                                space_sample_quantile, frametype, ctime,
                                calseq_sp_stare_max_t_duration,
                                n_calseq_sp_stare_fr, n_calseq_ca_stare_fr,
                                TIRS_T, 0.)
        i_sps += 1
        if i_sps == n_sps:
            break
        s0 = space_poss_slices[i_sps][0]

    # Case(s) after known HS:
    while n_HS > 0 and i_HS < n_HS and i_sps < n_sps:
        # If this is an immediately-pre-contact/maneuver HS, continue to the
        #  next HS, since there should be no cal sequences until after the
        #  next HS:
        if frametype[HS_slices[i_HS].stop+2] == TIRSp.CALTGT_POSS_POBS:
            i_HS += 1
            continue

        # Check for space aperture short scan (SASS):

        SASS_ib = HS_slices[i_HS].stop
        SASS_ct_b = ctime[SASS_ib-1]+TIRSp.ROIC_tau*0.5  # [s]
        SASS_ct_e = SASS_ct_b+35.94  # [s]
        SASS_ie = SASS_ib+np.searchsorted(ctime[SASS_ib:],
                                          SASS_ct_e)  # NumPy indexing

        if ctime[s0.start] > SASS_ct_b and ctime[s0.start] < SASS_ct_e:

            frametype[SASS_ib:SASS_ie] = TIRSp.SASS

            half_DNrange = np.mean(np.quantile(DN_c12_s4[SASS_ib:SASS_ie],
                                               [0.05,  0.9]))
            i1 = (np.nonzero(DN_c12_s4[SASS_ib:s0.start] > half_DNrange)[0][-1]+
                                                                         SASS_ib)
            i2 = (np.nonzero(DN_c12_s4[s0.stop:SASS_ie] > half_DNrange)[0][0]+
                                                                       s0.stop)
            sa_center_ct = (
                 sum([((ctime[i]-ctime[i-1])/(DN_c12_s4[i]-DN_c12_s4[i-1])*
                 (half_DNrange-DN_c12_s4[i])+ctime[i]) for i in [i1, i2]])*0.5)
            
            mirror_aostart_ctime.append(SASS_ct_b)  # [s]
            mirror_ang_offset.append((sa_center_ct-(SASS_ref_rel_ct+SASS_ct_b))*
                                     SASS_angspd)  # [deg]
           
            # Switch to the next cal-sequence space stare candidate, if any:
            i_sps += 1
            if i_sps == n_sps:
                break
            s0 = space_poss_slices[i_sps][0]
        else:
            pass  #### Add a non-zero default mirror angle here?

        # Info about the next HS, if any:
        if i_HS+1 < n_HS:
            i_of_HS_edge = HS_slices[i_HS+1].stop-1
        else:
            i_of_HS_edge = 999999999

        while s0.start <= i_of_HS_edge:
            n_fr = s0.stop-s0.start
            if (n_fr >= n_calseq_sp_stare_fr[0] and
                n_fr <= n_calseq_sp_stare_fr[1] and
                ctime[s0.stop-1]-ctime[s0.start] <
                        calseq_sp_stare_max_t_duration):
                # This cal-sequence space stare candidate has an acceptable
                #  frame length and time duration; vet it further:
                passed, frametype = vet_SAT1_sps(s0, space_poss_slices[i_sps],
                                DN_c12_s4, DN_c12_s6, deltaDN_calseq_info,
                                space_sample_quantile, frametype, ctime,
                                calseq_sp_stare_max_t_duration,
                                n_calseq_sp_stare_fr, n_calseq_ca_stare_fr,
                                TIRS_T, mirror_ang_offset[-1])
            i_sps += 1
            if i_sps == n_sps:
                break
            s0 = space_poss_slices[i_sps][0]

        i_HS += 1

    # Determine invalid or unusable apparent/possible space frames:
    ind_space0 = np.nonzero(np.abs(frametype) == TIRSp.SPACE_POSS)
    if len(ind_space0[0]) > 0:
        frametype, space_edge_t, space_block_ind = (
                    _find_invalid_calseqpart_fr_B(TIRSp.SPACE_POSS, frametype,
                                                  ind_space0, ctime, DN_c12_s4,
                                                  DN_c12_s6, use_min_DN=True))
        frametype[frametype == -TIRSp.SPACE_POSS] = TIRSp.SPACE_POSS
    else:
        space_edge_t, space_block_ind = (None, None)

    # Determine invalid or unusable apparent/possible caltgt frames:
    #  (must be called *after* the similar determination of space frames)
    ind_caltgt0 = np.nonzero(np.abs(frametype) == TIRSp.CALTGT_POSS)
    if len(ind_caltgt0[0]) > 0:
        frametype, caltgt_edge_t, caltgt_block_ind = (
                   _find_invalid_calseqpart_fr_B(TIRSp.CALTGT_POSS, frametype,
                                                 ind_caltgt0, ctime, DN_c12_s4,
                                                 DN_c12_s6))
        frametype[frametype == -TIRSp.CALTGT_POSS] = TIRSp.CALTGT_POSS
    else:
        caltgt_edge_t, caltgt_block_ind = (None, None)

    if len(mirror_aostart_ctime) == 0:
        mirror_aostart_ctime, mirror_ang_offset = (None, None)

    return (frametype, (space_edge_t, space_block_ind),
            (caltgt_edge_t, caltgt_block_ind),
            (mirror_aostart_ctime, mirror_ang_offset))


#--------------------------------------------------------------------------
def perform_categorization(sensor_IDval, bus_status_bfs, bus_events_bfs,
                             HS_slices, ctime, encoder_pos, obstgt_enc,
                             space_enc, caltgt_enc, DN_c12_s4, DN_c12_s6,
                             TIRS_T):
    """
       Identify calibration sequence portions (and skew frames)
         SAT1: without using motor encoder position values
         SAT2: primarily using motor encoder values
    """

#   Parameters:
#   -----------
#   sensor_IDval : int
#       Specifies which sensor is relevant, 1 = TIRS1, 2 = TIRS2
#   bus_status_bfs : ndarray (integer)
#       Array of values composed of bitflags, provides info related to bus status
#   bus_events_bfs : ndarray (integer)
#       Array of values composed of bitflags, provides info related to bus events
#   HS_slices : list
#       List of slices that indicate when the TIRS motor is powered off
#       ("hard stop", or HS)
#   ctime : ndarray (float64)
#       1D vector of continuous-time values (timestamps) [s]
#   encoder : ndarray (integer)
#       1D vector of motor encoder positions
#   DN_c12_s4 : ndarray (integer)
#       1D vector of ROIC DNs for channel 12, scene 4
#   DN_c12_s6 : ndarray (integer)
#       1D vector of ROIC DNs for channel 12, scene 6
#   obstgt_enc, space_enc, caltgt_enc : list of integers
#       Expected encoder positions for those frame types
#   DN_c12_s4 : ndarray (integer)
#       1D vector of ROIC DNs for channel 12, scene 4
#   DN_c12_s6 : ndarray (integer)
#       1D vector of ROIC DNs for channel 12, scene 6
#   TIRS_T : ndarray (float32)
#       TIRS thermistor temperatures

    # Initialize all frames as being nominal obstgt frame types:
    frametype = np.ma.empty((len(ctime),), dtype="int8")
    frametype[:] = TIRSp.OBSTGT_NOM

    oskew_b = None
    ctgtposs_b = None
    if sensor_IDval == 2:  # SAT2
        # Find any obvious-skew (oskew) frames (i.e., those that have a motor
        #  encoder position value outside of any of the expected frame types'
        #  ranges), then record their IDs in the 'frametype' array:
        oskew_enc = list(chain.from_iterable([obstgt_enc, caltgt_enc,
                                              space_enc]))
        some_oskew_fr, fr_ind = _find_via_encoder(encoder_pos, oskew_enc,
                                                  use_xor_neq=True)
        if some_oskew_fr:
            frametype[fr_ind] = TIRSp.SKEW_VIA_ENCODER
            oskew_b = ( frametype == TIRSp.SKEW_VIA_ENCODER )

        # Find any <possible> space frames, then record their IDs in the
        #  'frametype' array:
        some_possible_space_fr, fr_ind = _find_via_encoder(encoder_pos,
                                                           space_enc)
        if some_possible_space_fr:
            frametype[fr_ind] = TIRSp.SPACE_POSS

        # Find any <possible> caltgt frames, then record their IDs in the
        #  'frametype' array:
        some_possible_caltgt_fr, fr_ind = _find_via_encoder(encoder_pos,
                                                            caltgt_enc)
        if some_possible_caltgt_fr:
            frametype[fr_ind] = TIRSp.CALTGT_POSS
            ctgtposs_b = ( frametype == TIRSp.CALTGT_POSS )

    # Categorize all estimated payload-on-but-safed period frames as such (note
    #  that a minority of such frames may be recategorized later):
    during_PObS_period = get_bit_from_bitflags_v(bus_events_bfs, 5)
    frametype[during_PObS_period] = TIRSp.CALTGT_POSS_POBS

    # Set all motor-power-off frames as such, and categorize TIRS power-on
    #  "artifact" frames appropriately:
    for s in HS_slices:
        frametype[s] = TIRSp.HS_MPOFF

        ct_min = ctime[s.start]-46.  # 46 s (max time from TIRS power-on to
                                     #       reliably-"normal" DN values)
        i = s.start
        mcount = 0
        while ctime[i] >= ct_min and i > 0:
            i -= 1
            if np.abs(DN_c12_s4[i]-15000) > 10000:
                mcount += 1
        if mcount > 5:  # This seems to be a TIRS power-on case
            ib = i+1
            if ib == 1:  # Handle rare edge case
                ib = 0
            frametype[ib:s.start] = TIRSp.PON_ARTIFACTS

    # Determine space stare candidates from DN information:

    if sensor_IDval == 1:  # SAT1
        DN4_Tmoly_m, DN4_Tmoly_b = (-25.547, 23662.847)
        DN6_Tmoly_m, DN6_Tmoly_b = (-91.971, 45858.248)
        DN_Tmoly_nom_DNoffset = 200.
        i_Een_T_active = 1  # "telescope"
        i_Een_T_passive = 0  # "toroidal mirror"
    else:  # SAT2
        DN4_Tmoly_m, DN4_Tmoly_b = (-25.547, 23572.847)
        DN6_Tmoly_m, DN6_Tmoly_b = (-96.350, 45971.022)
        DN_Tmoly_nom_DNoffset = 200.
        i_Een_T_active = 2  # "filter module"
        i_Een_T_passive = 0  # "toroidal mirror"

    T_moly = TIRS_T[:,7]
    DN4_th0 = DN4_Tmoly_m*T_moly+DN4_Tmoly_b
    DN6_th0 = DN6_Tmoly_m*T_moly+DN6_Tmoly_b

      # "Representative" temperature difference, used to diagnose when the
      #  disruptive indirect-filter-heating radiometric event is occurring
      #  (near eclipse entrance). The associated DN offset below is added to
      #  the normal DN threshold in those cases.
    Een_T_ratio = TIRS_T[:,i_Een_T_active]/TIRS_T[:,i_Een_T_passive]
    Een_DN_Tmoly_nom_DNoffset = np.where(Een_T_ratio > 1.00081, 200., 0.)

      # At "very cold" operating temperatures seen only after an extended
      #  period of TIRS being powered off, the DN offset needs to be increased:
    cold_DN_Tmoly_nom_DNoffset = np.where(T_moly < 279., 200., 0.)

      # When operating temperatures are changing (usually increasing)
      #  rapidly enough (e.g., after a substantial period of TIRS being powered
      #  off, the DN offset needs to be increased (here, proportional to
      #  the 30-second change of T_moly):
    n_chk_s = 30
    ttmp = ((T_moly[n_chk_s:]-T_moly[:-n_chk_s])/
            (ctime[n_chk_s:]-ctime[:-n_chk_s]))
    Tdiff_chk = np.insert(ttmp, 0, [0. for x in range(n_chk_s)])
    ttmp = np.abs(Tdiff_chk)
    msk = ( ttmp > 1.e-3 )
    cold_DN_Tmoly_nom_DNoffset[msk] = 8.e4*ttmp[msk]

    DN_threshold = (DN4_th0+DN_Tmoly_nom_DNoffset+Een_DN_Tmoly_nom_DNoffset+
                    cold_DN_Tmoly_nom_DNoffset)
    space_s4_tnom_b = ( DN_c12_s4 < DN_threshold )

    DN_threshold = (DN6_th0+DN_Tmoly_nom_DNoffset+Een_DN_Tmoly_nom_DNoffset+
                    cold_DN_Tmoly_nom_DNoffset)
    space_s6_tnom_b = ( DN_c12_s6 < DN_threshold )

      # Both scenes (or equivalently in this case, ROICs) need to satisfy
      #  the criterion -- in order to help remain robust against internally- or
      #  externally-induced electrical noise:
    space_tnom_frames_b = np.logical_and(space_s4_tnom_b, space_s6_tnom_b)

    space_poss_frames_b = space_tnom_frames_b
    if oskew_b is not None:  # Mask any obviously-skew frames (encoder-based)
        space_poss_frames_b[oskew_b] = False
    if ctgtposs_b is not None:  # Mask any possible caltgt frames (encoder-based)
        space_poss_frames_b[ctgtposs_b] = False
    labels, n_space_tmp = scipy.ndimage.label(space_poss_frames_b)
    space_poss_slices = scipy.ndimage.find_objects(labels)

    if sensor_IDval == 1:  # SAT1
        frametype, space_result_t, caltgt_result_t, mirror_t = sift_SAT1_data(
                  frametype, space_poss_slices, bus_status_bfs, bus_events_bfs,
                  HS_slices, ctime, encoder_pos, DN_c12_s4, DN_c12_s6, TIRS_T)
    else:  # SAT2
        frametype, space_result_t, caltgt_result_t, mirror_t = sift_SAT2_data(
                  frametype, space_poss_slices, bus_status_bfs, bus_events_bfs,
                  HS_slices, ctime, encoder_pos, DN_c12_s4, DN_c12_s6, TIRS_T)

    # Remove any single OBSTGT_NOM frames sandwiched between other frame types:
    tmp0 = (frametype == TIRSp.OBSTGT_NOM)
    tmp = [i for i in range(1,len(tmp0)-1)
                                if list(tmp0[i-1:i+2]) == [False, True, False]]
    if len(tmp) > 0:
        frametype[tmp] = TIRSp.SKEW_VIA_DISP

    return (frametype, space_result_t, caltgt_result_t, mirror_t)


#--------------------------------------------------------------------------
def read_prepped_payload_tlm(ds_l, n_atrack_total, n_atrack_tl):
    """Read prepped payload telemetry file data."""

    # Read in data, concatenating arrays along the "atrack" dimension:
    ip = 0
    tmpdat = {}
    for i_ds in range(len(ds_l)):
        ds = ds_l[i_ds]
        if ip == 0:
            is_tlm = ( "dummy_var_for_dims" not in
                                          ds.groups["TIRS_L0_Data"].variables )
            if is_tlm:
                vars_to_read = ["ctime", "encoder_pos", "PDU_12V_I",
                                "adj_srd_DN", "TIRS_T"]
                vars_to_autoproc = ["encoder_pos", "PDU_12V_I", "TIRS_T"]
            else:
                vars_to_read = ["ctime"]

            for vn in vars_to_read:
                v_dtype = ds.groups["TIRS_L0_Data"].variables[vn].dtype
                v_shape = ds.groups["TIRS_L0_Data"].variables[vn].shape
                if len(v_shape) == 1:  # (atrack)
                    concat_shape = (n_atrack_total,)
                elif len(v_shape) == 2:  # (atrack, m)
                    concat_shape = (n_atrack_total, v_shape[1])
                elif len(v_shape) == 3:  # (m, n, atrack)
                    concat_shape = (v_shape[0], v_shape[1], n_atrack_total)
                tmpdat[vn] = np.ma.empty(concat_shape, dtype=v_dtype)

        this_n_atrack, ia_b, ia_e = n_atrack_tl[i_ds]
        tmpdat["ctime"][ip:ip+this_n_atrack] = (
                        ds.groups["TIRS_L0_Data"].variables["ctime"][ia_b:ia_e])

        if is_tlm:
            for vn in vars_to_autoproc:
                tmpdat[vn][ip:ip+this_n_atrack] = (
                             ds.groups["TIRS_L0_Data"].variables[vn][ia_b:ia_e])
            tmpdat["adj_srd_DN"][:,:,ip:ip+this_n_atrack] = (
                                      ds.groups["TIRS_L0_Data"].
                                         variables["adj_srd_DN"][:,:,ia_b:ia_e])
        ip += this_n_atrack

    # Sort arrays by ctime:
    odat = {}
    sorted_inds = np.argsort(tmpdat["ctime"], kind="stable")
    odat["ctime"] = tmpdat["ctime"][sorted_inds]
    if is_tlm:
        for vn in vars_to_autoproc:
            odat[vn] = tmpdat[vn][sorted_inds]

    # Chose channel 12 (here, Python index 12 of [0-63]) because
    #   * it is in the "deepest" portion of the mid-IR atmospheric window
    #     which should maximize contrast with the cold space frames
    #   * it is good for scenes 1, 4, 5, 6, and 8 (for both TIRS)
    # Chose scenes 4 and 6 (here, Python indices 3, 5 of [0-7]) because
    #   * they are good for both TIRS
    #   * they are physically separated to reduce the chances that an energetic
    #     particle hit will affect both over the same short time interval
    if is_tlm:
        odat["adj_srd_DN_c12_s4"] = tmpdat["adj_srd_DN"][12,3,sorted_inds]
        odat["adj_srd_DN_c12_s6"] = tmpdat["adj_srd_DN"][12,5,sorted_inds]

    return odat


#--------------------------------------------------------------------------
def read_prepped_bus_tlm(ds_l, n_atrack_total, n_atrack_tl):
    """Read prepped bus telemetry file data."""

    # Read in data, concatenating arrays along the "atrack" dimension:
    ip = 0
    tmpdat = {}
    for i_ds in range(len(ds_l)):
        ds = ds_l[i_ds]
        if ip == 0:
            for vn in bus_vars_to_read:
                v_dtype = ds.groups["SCbus_L0_Data"].variables[vn].dtype
                v_shape = ds.groups["SCbus_L0_Data"].variables[vn].shape
                if len(v_shape) == 1:
                    concat_shape = (n_atrack_total,)
                elif len(v_shape) == 3:
                    concat_shape = (v_shape[0], v_shape[1], n_atrack_total)
                tmpdat[vn] = np.ma.empty(concat_shape, dtype=v_dtype)

        this_n_atrack, ia_b, ia_e = n_atrack_tl[i_ds]
        for vn in bus_vars_to_read:
            tmpdat[vn][ip:ip+this_n_atrack] = (
                            ds.groups["SCbus_L0_Data"].variables[vn][ia_b:ia_e])
        ip += this_n_atrack

    odat = {}
    sorted_inds = np.argsort(tmpdat["ctime"], kind="stable")
    for vn in bus_vars_to_read:
        odat[vn] = tmpdat[vn][sorted_inds]

    return odat


#--------------------------------------------------------------------------
def read_prepped_orbit_reconst(ds_l, n_atrack_total, n_atrack_tl):
    """Read prepped orbit reconstruction file data."""

    # Read in data, concatenating arrays along the "atrack" dimension:
    ip = 0
    tmpdat = {}
    for i_ds in range(len(ds_l)):
        ds = ds_l[i_ds]
        if ip == 0:
            for vn in orbit_vars_to_read:
                v_dtype = ds.groups["Orbit_L0_Data"].variables[vn].dtype
                v_shape = ds.groups["Orbit_L0_Data"].variables[vn].shape
                if len(v_shape) == 1:
                    concat_shape = (n_atrack_total,)
                elif len(v_shape) == 3:
                    concat_shape = (v_shape[0], v_shape[1], n_atrack_total)
                tmpdat[vn] = np.ma.empty(concat_shape, dtype=v_dtype)

        this_n_atrack, ia_b, ia_e = n_atrack_tl[i_ds]
        for vn in orbit_vars_to_read:
            tmpdat[vn][ip:ip+this_n_atrack] = (
                            ds.groups["Orbit_L0_Data"].variables[vn][ia_b:ia_e])
        ip += this_n_atrack

    odat = {}
    sorted_inds = np.argsort(tmpdat["ctime"], kind="stable")
    for vn in orbit_vars_to_read:
        odat[vn] = tmpdat[vn][sorted_inds]

    return odat


#--------------------------------------------------------------------------
def _find_any_usable_calseqs(ctime, space_t, caltgt_t, frametype):
    max_n_cspart_elems = 24
    w_frametype = frametype.copy()

    space_edge_t, space_block_ind = space_t
    caltgt_edge_t, caltgt_block_ind = caltgt_t

    if (space_edge_t is None) or (caltgt_edge_t is None):
        return (frametype, -9999)  # Cannot refine more (no cal-seqs)
    
    # Determine the correspondence of blocks to valid calibration sequences.
    #  A valid calibration sequence is one where the leading caltgt-block edge
    #  is < 12 seconds after a trailing space-block edge. 
    ia_space_ld_edge, ia_space_tr_edge = space_edge_t
    ia_caltgt_ld_edge, ia_caltgt_tr_edge = caltgt_edge_t
    ts_space_tr_e = ctime[ia_space_tr_edge]
    ts_caltgt_ld_e = ctime[ia_caltgt_ld_edge]
    calseq_ind_tl = []
    for i_caltgtb in range(len(ia_caltgt_ld_edge)):
         testv = ts_caltgt_ld_e[i_caltgtb]-ts_space_tr_e
         s_ind = np.nonzero((testv > 0.) & (testv < 12.))
         if len(s_ind[0]) > 0:
             calseq_ind_tl.append((s_ind[0][0], i_caltgtb))
    n_usable_calseqs = len(calseq_ind_tl)
    if n_usable_calseqs == 0:
        return (w_frametype, -9999)  # Cannot refine more (no cal-seqs)

    for i_spaceb, i_caltgtb in calseq_ind_tl:
        w_frametype[space_block_ind[i_spaceb]] = TIRSp.SPACE_NOM
        w_frametype[caltgt_block_ind[i_caltgtb]] = TIRSp.CALTGT_NOM

    return (w_frametype, n_usable_calseqs)


#--------------------------------------------------------------------------
def apply_dat_edges_by_ct_to_scidat_bool(dat, scidat, dat_ct_margin_t,
                            dat_inds=None, dat_slices=None, bool_v=None):
    """
       Given a 1-D list or NumPy array ('dat_inds' -OR- 'dat_slices'; non-payload
       data stream indices or slices, respectively), non-payload and
       payload data stream dictionaries ('dat' and 'scidat', respectively),
       and ctime margins to use before and after dat ctime edges, return a
       sequence of boolean values (valid at each scidat ctime) indicating the
       scidat stream frames that the dat indices/slices correspond to in
       ctime space.  Optionally, 'bool_v' specifies initial values for the
       output.
    """

    if dat_inds is None and dat_slices is None:
        raise RuntimeError("must provide non-payload data info; "
                           "indices -OR- slices")
    elif dat_inds is not None and dat_slices is not None:
        raise RuntimeError("must not provide BOTH non-payload data indices "
                           "and slices")

    if dat_inds is not None:
        if len(dat_inds) == 0:
            if bool_v is None:
                return dat_inds  # Return empty input
            else:
                return bool_v.copy()  # Return copy of "output initial values"

        if dat_inds[0] == 0:
            dat_ct_tl = [(dat["ctime"][0]-1.,#@data_ct_margin_t[0],
                             dat["ctime"][0])]
            tmp_tl = [(dat["ctime"][i-1]+dat_ct_margin_t[0],
                       dat["ctime"][i]+dat_ct_margin_t[1]) for i in
                       dat_inds[1:]]
            dat_ct_tl = dat_ct_tl+tmp_tl
        else:
            dat_ct_tl = [(dat["ctime"][i-1]+dat_ct_margin_t[0],
                             dat["ctime"][i]+dat_ct_margin_t[1]) for i in
                             dat_inds]

    if dat_slices is not None:
        if len(dat_slices) == 0:
            if bool_v is None:
                return []  # Return empty input
            else:
                return bool_v.copy()  # Return copy of "output initial values"

        dat_ct_tl = [(dat["ctime"][s.start]+dat_ct_margin_t[0],
                      dat["ctime"][s.stop-1]+dat_ct_margin_t[1]) for s in
                         dat_slices]

    if bool_v is None:
        bool_tmp = np.full(scidat["ctime"].shape, False, dtype="bool")
    else:
        bool_tmp = bool_v.copy()
    for ct_b, ct_e in dat_ct_tl:
        bool_tmp[(scidat["ctime"] >= ct_b) & (scidat["ctime"] <= ct_e)] = True
    return bool_tmp


#--------------------------------------------------------------------------
def apply_slices_to_bit_in_bitflags_v(dat, scidat, dat_ct_margin_t_inc,
                                      slice_msk, offset, int_v):
    lbls, _ = scipy.ndimage.label(slice_msk)
    slices = [s[0] for s in scipy.ndimage.find_objects(lbls)]
    bool_tmp = apply_dat_edges_by_ct_to_scidat_bool(dat, scidat,
                                   dat_ct_margin_t_inc, dat_slices=slices)
    if len(bool_tmp) > 0:
        return apply_bit_to_bitflags_v(offset, bool_tmp, int_v)
    else:
        return int_v.copy()


#--------------------------------------------------------------------------
def categorize_L0_payload_tlm(leap_s_info, ext_cfg_d=None):
    """Driver routine."""

    cfg_d = {}
    if ext_cfg_d is None:
        cfg_d["ref_L0_payload_paths"] = os.environ["L0_PAYLOAD_FPATHS"]
        if "PLD_TMPCUR_FPATH" in os.environ:
            cfg_d["L0_payload_paths"] = os.environ["PLD_TMPCUR_FPATH"]
        else:
            cfg_d["L0_payload_paths"] = os.environ["L0_PAYLOAD_FPATHS"]

        cfg_d["ref_L0_bus_paths"] = os.environ["L0_BUS_FPATHS"]
        if "BUS_TMPCUR_FPATH" in os.environ:
            cfg_d["L0_bus_paths"] = os.environ["BUS_TMPCUR_FPATH"]
        else:
            cfg_d["L0_bus_paths"] = os.environ["L0_BUS_FPATHS"]

        cfg_d["output_filespecs_fpath"] = os.path.join(
                                              os.environ["ANCILLARY_DATA_DIR"],
                                         "categorization_L0_pld_filespecs.json")
        cfg_d["output_dir"] = os.environ["OUTPUT_DIR"]
        cfg_d["product_fullver"] = os.environ["PRODUCT_FULLVER"]

        tmp_l = os.environ["GRANULE_START_ID_END"].split('_')
        cfg_d["ctime_range_to_proc"] = (float(tmp_l[0])-16.99*hour_to_s,
                                        float(tmp_l[2])+1.59*hour_to_s)
    else:
        cfg_d = ext_cfg_d

    # Get some top-level info about the input file(s):
    input_b_fpaths = [x.strip() for x in cfg_d["L0_bus_paths"].split("||")]
    ref_input_b_fpaths = [x.strip() for x in
                          cfg_d["ref_L0_bus_paths"].split("||")]
    b_info = general_tlm_file_info(input_b_fpaths, "bus",
                                   start_ctime=cfg_d["ctime_range_to_proc"][0],
                                   end_ctime=cfg_d["ctime_range_to_proc"][1])

    input_r_fpaths = [x.strip() for x in cfg_d["L0_orbit_paths"].split("||")]
    ref_input_r_fpaths = [x.strip() for x in
                          cfg_d["ref_L0_orbit_paths"].split("||")]
    r_info = general_tlm_file_info(input_r_fpaths, "orbit_reconst",
                                   start_ctime=cfg_d["ctime_range_to_proc"][0],
                                   end_ctime=cfg_d["ctime_range_to_proc"][1])

    input_p_fpaths = [x.strip() for x in cfg_d["L0_payload_paths"].split("||")]
    ref_input_p_fpaths = [x.strip() for x in
                          cfg_d["ref_L0_payload_paths"].split("||")]
    p_info = general_tlm_file_info(input_p_fpaths, "payload",
                                   start_ctime=cfg_d["ctime_range_to_proc"][0],
                                   end_ctime=cfg_d["ctime_range_to_proc"][1])

    # Read relevant fields from the input file(s); then close the files:
    scidat = read_prepped_payload_tlm(p_info["ds_l"], p_info["n_atrack_total"],
                                      p_info["n_atrack_tl"])
    for ds in p_info["ds_l"]:
        ds.close()

    busdat = read_prepped_bus_tlm(b_info["ds_l"], b_info["n_atrack_total"],
                                  b_info["n_atrack_tl"])
    for ds in b_info["ds_l"]:
        ds.close()

    orbdat = read_prepped_orbit_reconst(r_info["ds_l"], r_info["n_atrack_total"],
                                        r_info["n_atrack_tl"])
    for ds in r_info["ds_l"]:
        ds.close()

    o_data = {}
    bus_status_bfs = np.ma.zeros((p_info["n_atrack_total"],), dtype="uint8")
    bus_events_bfs = np.ma.zeros((p_info["n_atrack_total"],), dtype="uint16")
    frametype = np.ma.ones((p_info["n_atrack_total"],), dtype="int8")*(-99)

    busdat_ct_margin_t_exc = (0.1, 0.1)  # [s] (begin+x, end-x); exclusive
                                         #  best if < minimum bus tlm cadence?
    busdat_ct_margin_t_inc = (0.1, 0.1)  # [s] (begin-x, end+x); inclusive
                                         #  best if < minimum bus tlm cadence?
    orbdat_ct_margin_t_inc = (0.1, 0.1)  # [s] (begin-x, end+x); inclusive
                                         #  best if < minimum orbrec cadence?
    scidat_ct_margin_t_inc = (0.1, 0.1)  # [s] (begin-x, end+x); inclusive
                                         #  best if < minimum bus tlm cadence?

   #--- Set bus status bitflags field:

    # Record bus telemetry gaps (b3):
#^test    busdat_gap_threshold_t = (6., 5.4)  # [s]
    busdat_gap_threshold_t = (0.19, 5.4)  # [s]
    busdat_gap_dcheck = np.append(np.array(busdat_gap_threshold_t[0]),
                                  np.diff(busdat["ctime"]))
    busdat_gap_inds = np.nonzero(
                              busdat_gap_dcheck > busdat_gap_threshold_t[1])[0]
    bool_tmp = apply_dat_edges_by_ct_to_scidat_bool(busdat, scidat,
                               busdat_ct_margin_t_exc, dat_inds=busdat_gap_inds)
    if len(bool_tmp) > 0:
        bus_status_bfs = apply_bit_to_bitflags_v(3, bool_tmp, bus_status_bfs)

    # Record when REFS is invalid (b2):
    slice_msk = ( busdat["REFS_is_valid"].data == 0 )
    bus_status_bfs = apply_slices_to_bit_in_bitflags_v(busdat, scidat,
                          busdat_ct_margin_t_inc, slice_msk, 2, bus_status_bfs)

    # Record when ATT_DET is invalid (b1):
    slice_msk = ( busdat["ATT_DET_is_valid"].data == 0 )
    bus_status_bfs = apply_slices_to_bit_in_bitflags_v(busdat, scidat,
                          busdat_ct_margin_t_inc, slice_msk, 1, bus_status_bfs)

    # Record when REFS_ESM is invalid (b0):
    slice_msk = ( busdat["REFS_ESM_is_valid"].data == 0 )
    bus_status_bfs = apply_slices_to_bit_in_bitflags_v(busdat, scidat,
                          busdat_ct_margin_t_inc, slice_msk, 0, bus_status_bfs)

    # Record when payload is powered off (b4):
    slice_msk = ( busdat["payload_is_on"].data == 0 )
    bus_status_bfs = apply_slices_to_bit_in_bitflags_v(busdat, scidat,
                          busdat_ct_margin_t_inc, slice_msk, 4, bus_status_bfs)

   #--- Set bus events bitflags field:

    # Record when pass is ascending (b0):
    slice_msk = ( orbdat["velocity_wrt_ECI_3"].data >= 0. )
    bus_events_bfs = apply_slices_to_bit_in_bitflags_v(orbdat, scidat,
                          orbdat_ct_margin_t_inc, slice_msk, 0, bus_events_bfs)

    # Record when in moon umbra (b1):
    slice_msk = ( busdat["SC_in_moon_umbra"].data == 1 )
    bus_events_bfs = apply_slices_to_bit_in_bitflags_v(busdat, scidat,
                          busdat_ct_margin_t_inc, slice_msk, 1, bus_events_bfs)

    # Record when in Earth umbra (b2):
    slice_msk = ( orbdat["SC_in_Earth_umbra"].data == 1 )
    bus_events_bfs = apply_slices_to_bit_in_bitflags_v(orbdat, scidat,
                          orbdat_ct_margin_t_inc, slice_msk, 2, bus_events_bfs)

    # Record when radio is tx (b3):
    slice_msk = ( busdat["radio_is_tx"].data == 1 )
    bus_events_bfs = apply_slices_to_bit_in_bitflags_v(busdat, scidat,
                          busdat_ct_margin_t_inc, slice_msk, 3, bus_events_bfs)

    # Record when radio is active (tx and/or rx) (b4):
    radio_on_slice_msk = ( (busdat["radio_is_tx"].data == 1) |
                           (busdat["radio_rx_lock"].data == 1) )
    bus_events_bfs = apply_slices_to_bit_in_bitflags_v(busdat, scidat,
                 busdat_ct_margin_t_inc, radio_on_slice_msk, 4, bus_events_bfs)

    # Record when attitude differs from commanded attitude (b6+):

    diff_q1 = busdat["q_scbody_wrt_ECI_1"]-busdat["cmd_q_scbody_wrt_ECI_1"]
    diff_q2 = busdat["q_scbody_wrt_ECI_2"]-busdat["cmd_q_scbody_wrt_ECI_2"]
    diff_q3 = busdat["q_scbody_wrt_ECI_3"]-busdat["cmd_q_scbody_wrt_ECI_3"]
    diff_q4 = busdat["q_scbody_wrt_ECI_4"]-busdat["cmd_q_scbody_wrt_ECI_4"]
    slice_msk = ( (np.abs(diff_q1) > 0.015) | (np.abs(diff_q2) > 0.015) |
                  (np.abs(diff_q3) > 0.015) | (np.abs(diff_q4) > 0.015) )
    lbls, _ = scipy.ndimage.label(slice_msk)
    att_diff_slices = [s[0] for s in scipy.ndimage.find_objects(lbls)]

       # Only consider primary-commanded-vector changes over a short-enough
       #  time window (this is a mitigation to deal with the presence of longer
       #  bus telemetry gaps):
    pri_cmd_ref_diff = np.diff(busdat["pri_cmd_ref_dir"])
    pri_cmd_ref_chg_t = [(i+1, tuple(busdat["pri_cmd_ref_dir"].data[i:i+2]))
                         for i in np.nonzero(pri_cmd_ref_diff)[0] if
                         busdat["ctime"][i+1]-busdat["ctime"][i] < 10.5]

       # Categorize attitude deviation slices:
    att_diff_slice_cat_d = {0: [], 1: [], 2: [], 3: [], 4: [], 5: []}
    prior_cat = 0
    for s in att_diff_slices:
        cat = 0  # Default (unknown category)
        for ii, chg_t in pri_cmd_ref_chg_t:
            if ii >= s.start and ii <= s.stop-1:
                cat = pri_cmd_chg_t_ref_d[chg_t]
        if np.any(busdat["pri_cmd_vec_wrt_SBF_1"][s] < -0.95):
            cat = 4  # point -X (GPS receiver) to zenith
        elif ((prior_cat == 4) and
              np.any(busdat["pri_cmd_vec_wrt_SBF_1"][s] > -0.05)):  #@ > 0.95
            cat = 2  # slew back to LVLH (from cat = 4)
        elif np.any(busdat["MSA_action_occurring"][s] == 1):
            cat = 5  # MSA action is occcurring
        att_diff_slice_cat_d[cat].append(s)
        prior_cat = cat
       # Record attitude deviation slices (b6+):
    for key in att_diff_slice_cat_d:
        slices = att_diff_slice_cat_d[key]
        if len(slices) > 0:
            bool_tmp = apply_dat_edges_by_ct_to_scidat_bool(busdat, scidat,
                                      busdat_ct_margin_t_inc, dat_slices=slices)
            if len(bool_tmp) > 0:
                bus_events_bfs = apply_bit_to_bitflags_v(6+key, bool_tmp,
                                                         bus_events_bfs)

    # Record when within a payload-on-but-safed bus maneuver period (b5; e.g.,
    #  contacts, GPS-zenith-pointing) -- approximately equivalent to:
    #      start of BusMacro44 motor hard stop (HS) to start of BusMacro45 HS
    #
    # Contact period sequence of events:
    #    {Macro44 (HS + begin cal stare; 90 s before scheduled contact),
    #     wait,
    #     bus slew (begins 60 s before scheduled contact),
    #     wait,
    #     radio tx/rx,
    #     bus slew (begins very near the scheduled contact end time),
    #     wait,
    #     Macro45 (HS + run MT 3, then MT 6; 75 s after scheduled contact end)}
    #
    # GPS-zenith-pointing sequence of events:
    #    {Macro44 (HS + begin cal stare; 15 s before scheduled maneuver),
    #     wait,
    #     bus slew (begins at scheduled maneuver start time),
    #     wait,
    #     bus slew (begins at scheduled maneuver end time),
    #     wait,
    #     Macro45 (HS + run MT 3, then MT 6; 60 s after scheduled maneuver end)}
    #
       # First, estimate contact period(s) using bus attitude events:
    PObS_busmnv_ct_bounds = []
    for slew_to_contact in att_diff_slice_cat_d[1]:
        ct_beg = busdat["ctime"][slew_to_contact.start]  # [s]
        for slew_to_LVLH in att_diff_slice_cat_d[2]:
            ct_end = busdat["ctime"][slew_to_LVLH.start]  # [s]
            ct_delta = ct_end-ct_beg  # [s]
            if ct_delta < 0.:
                continue
            elif ct_delta <= 660.:  # 10-minute max contact + 60 s
                PObS_busmnv_ct_bounds.append((ct_beg-30., ct_end+75.))
            else:
                break
        # Second, estimate GPS-zenith-pointing period(s) using bus att events:
    for slew_to_zenith in att_diff_slice_cat_d[4]:
        ct_beg = busdat["ctime"][slew_to_zenith.start]  # [s]
        for slew_to_LVLH in att_diff_slice_cat_d[2]:
            ct_end = busdat["ctime"][slew_to_LVLH.stop-1]  # [s]
            ct_delta = ct_end-ct_beg  # [s]
            if ct_delta < 0.:
                continue
            elif ct_delta > 720. and ct_delta <= 1500.:  # 12 to 25 minutes
#@                PObS_busmnv_ct_bounds.append((ct_beg-15., ct_end+60.))
                PObS_busmnv_ct_bounds.append((ct_beg-15., ct_end+5.))
            else:
                break

       # Next try to determine contact period(s) using motor hard stops:
    slice_msk = ( scidat["PDU_12V_I"] < 0.008 )
    lbls, _ = scipy.ndimage.label(slice_msk)
    HS_slices = [s[0] for s in scipy.ndimage.find_objects(lbls)]
    HS_sl_ct_t = [(scidat["ctime"][sl.start], scidat["ctime"][sl.stop-1]) for
                                                               sl in HS_slices]
    PObS_busmnv_ct_bounds_HS = []
    if len(HS_slices) > 1:
        for isl in range(len(HS_slices)-1):
            ib_A, ib_B = (HS_slices[isl].start, HS_slices[isl+1].start)
            radio_is_tx = get_bit_from_bitflags_v(bus_events_bfs[ib_A+1:ib_B],
                                                  3)
            bus_slew_to_catt = get_bit_from_bitflags_v(
                                                bus_events_bfs[ib_A+1:ib_B], 7)
            bus_slew_to_nxzen = get_bit_from_bitflags_v(
                                               bus_events_bfs[ib_A+1:ib_B], 10)
            tmp = scidat["ctime"][ib_B]-scidat["ctime"][ib_A]
              # (90 s) + (3-minute min contact duration) + (75 s) = 345 s
              # (90 s) + (10-minute max contact duration) + (75 s) = 765 s
            HS_intv_okay_for_contact = ( (tmp >= 345.) & (tmp <= 765.) )
              # (15 s) + (12-minute min GPS-zenith-point) + (60 s) = 795 s
              # (15 s) + (25-minute max GPS-zenith-point) + (60 s) = 1575 s
            HS_intv_okay_for_nxzen = ( (tmp >= 795.) & (tmp <= 1575.) )

            if ((np.any(radio_is_tx) or np.any(bus_slew_to_catt)) and
                                                     HS_intv_okay_for_contact):
                PObS_busmnv_ct_bounds_HS.append((scidat["ctime"][ib_A],
                                               scidat["ctime"][ib_B]))
            elif (np.any(bus_slew_to_nxzen) and HS_intv_okay_for_nxzen):
                PObS_busmnv_ct_bounds_HS.append((scidat["ctime"][ib_A],
                                               scidat["ctime"][ib_B]))

    if len(PObS_busmnv_ct_bounds) > 0 and len(PObS_busmnv_ct_bounds_HS) > 0:
        tmp_ctb = []
        tmp_HSctb_used = [False for x in PObS_busmnv_ct_bounds_HS]
        for b_ctb_t in PObS_busmnv_ct_bounds:
            ctb_0, ctb_1 = b_ctb_t  # Default
            for ih, h_ctb_t in enumerate(PObS_busmnv_ct_bounds_HS):
                if np.abs(b_ctb_t[0]-h_ctb_t[0]) < 10.:
                    ctb_0 = min(b_ctb_t[0], h_ctb_t[0])
                    tmp_HSctb_used[ih] = True
                if np.abs(b_ctb_t[1]-h_ctb_t[1]) < 10.:
                    ctb_1 = max(b_ctb_t[1], h_ctb_t[1])
                    tmp_HSctb_used[ih] = True
            tmp_ctb.append((ctb_0, ctb_1))
           # Add any remaining HS-indicated event(s):
        for b, h_ctb_t in zip(tmp_HSctb_used, PObS_busmnv_ct_bounds_HS):
            if not b:
                tmp_ctb.append(h_ctb_t)
        PObS_busmnv_ct_bounds = tmp_ctb
    elif len(PObS_busmnv_ct_bounds) > 0:
        pass  # Already have needed info in 'PObS_busmnv_ct_bounds'
    elif len(PObS_busmnv_ct_bounds_HS) > 0:
        PObS_busmnv_ct_bounds = [x for x in PObS_busmnv_ct_bounds_HS]
    else:
        pass  # Nothing can be done?

    if len(PObS_busmnv_ct_bounds) > 0:
        bool_tmp = np.full(scidat["ctime"].shape, False, dtype="bool")
        for ct_b, ct_e in PObS_busmnv_ct_bounds:
            bool_tmp[(scidat["ctime"] >= ct_b-scidat_ct_margin_t_inc[0]) &
                     (scidat["ctime"] <= ct_e+scidat_ct_margin_t_inc[1])] = True
        bus_events_bfs = apply_bit_to_bitflags_v(5, bool_tmp, bus_events_bfs)

   #--- Determine type of each frame:

    is_tlm = ( "encoder_pos" in scidat )  # versus prototlm

    if is_tlm:
        if p_info["sensor_IDval"] == 1:
            ep_min, ep_max = TIRSp.obstgt_enc_pos_range[p_info["sensor_IDidx"]]
            obstgt_enc_pos_l = list(range(ep_min, ep_max+1))
            ep_min, ep_max = TIRSp.space_enc_pos_range[p_info["sensor_IDidx"]]
            space_enc_pos_l = list(range(ep_min, ep_max+1))
            ep_min, ep_max = TIRSp.caltgt_enc_pos_range[p_info["sensor_IDidx"]]
            caltgt_enc_pos_l = list(range(ep_min, ep_max+1))

            frametype, space_t, caltgt_t, mirror_t = (
                   perform_categorization(
                    p_info["sensor_IDval"], bus_status_bfs, bus_events_bfs,
                    HS_slices, scidat["ctime"], scidat["encoder_pos"],
                    obstgt_enc_pos_l, space_enc_pos_l, caltgt_enc_pos_l,
                    scidat["adj_srd_DN_c12_s4"], scidat["adj_srd_DN_c12_s6"],
                    scidat["TIRS_T"]))

        elif p_info["sensor_IDval"] == 2:
            ep_min, ep_max = TIRSp.obstgt_enc_pos_range[p_info["sensor_IDidx"]]
            obstgt_enc_pos_l = list(range(ep_min, ep_max+1))
            ep_min, ep_max = TIRSp.space_enc_pos_range[p_info["sensor_IDidx"]]
            space_enc_pos_l = list(range(ep_min, ep_max+1))
            ep_min, ep_max = TIRSp.caltgt_enc_pos_range[p_info["sensor_IDidx"]]
            caltgt_enc_pos_l = list(range(ep_min, ep_max+1))

            frametype, space_t, caltgt_t, mirror_t = (
                   perform_categorization(
                    p_info["sensor_IDval"], bus_status_bfs, bus_events_bfs,
                    HS_slices, scidat["ctime"], scidat["encoder_pos"],
                    obstgt_enc_pos_l, space_enc_pos_l, caltgt_enc_pos_l,
                    scidat["adj_srd_DN_c12_s4"], scidat["adj_srd_DN_c12_s6"],
                    scidat["TIRS_T"]))

        else:
            raise ValueError(
                f"sensor ID value of {p_info['sensor_IDval']} is not supported")
  
        # Determine any usable calibration sequences, finalize 'frametype'
        #  array:
        n_usable_calseqs = -9999
        if np.any(frametype):
            frametype, n_usable_calseqs = _find_any_usable_calseqs(
                scidat["ctime"], space_t, caltgt_t, frametype)

    else:
        n_usable_calseqs = -9999

    o_data["ctime"] = scidat["ctime"]
    o_data["scipkt_frame_type"] = frametype
    o_data["scipkt_bus_status_bitflags"] = bus_status_bfs
    o_data["scipkt_bus_events_bitflags"] = bus_events_bfs

    if mirror_t[0] is None: 
        # Need to have at least one value (default of zero), even if no mirror
        #  info is found:
        o_data["mirror_aostart_ctime"] = np.array((scidat["ctime"][0],))
        o_data["mirror_ang_offset"] = np.array((0.,))
    else:
        o_data["mirror_aostart_ctime"] = np.array(mirror_t[0])
        o_data["mirror_ang_offset"] = np.array(mirror_t[1])

    #== Determine/set some global/group attribute values:

    product_full_version = cfg_d["product_fullver"].strip()

    ctime_coverage = np.array([np.amin(scidat["ctime"]),
                               np.amax(scidat["ctime"])])  # [s]
    UTC_DT, _ = ctime_to_UTC_DT(ctime_coverage, 's', leap_s_info)
    UTC_coverage = ([UTC_DT[0].strftime("%Y-%m-%dT%H:%M:%S.%f"),
                     UTC_DT[1].strftime("%Y-%m-%dT%H:%M:%S.%f")])

    with open(L1_fpaths.scipkg_prdgitv_fpaths[3], 'r') as in_f:
        line_parts = in_f.readline().split('(', maxsplit=1)
        this_pkg_provenance = "{}{} ( {}".format(line_parts[0],
                                                 product_full_version,
                                                 line_parts[1].strip())

    with open(L1_fpaths.scipkg_version_fpath, 'r') as in_f:
        proc_algID = in_f.readline().strip()

    ipf_l = ([os.path.basename(x) for x in ref_input_p_fpaths]+
             [os.path.basename(x) for x in ref_input_b_fpaths]+
             [os.path.basename(x) for x in ref_input_r_fpaths])
    input_product_files = ", ".join(ipf_l)

      # Finish output filepaths:
    outp_fn_body = f"prefire_0{p_info['sensor_IDval']:1d}_pld_cat"
    now_UTC_DT = datetime.datetime.now(datetime.timezone.utc)
    now_UTC_strrep = now_UTC_DT.strftime("%Y-%m-%dT%H:%M:%S.%f")
    outp_fn_suffix = \
         f"_{UTC_DT[0]:%Y%m%d%H%M%S}_{UTC_DT[1]:%Y%m%d%H%M%S}_{now_UTC_DT:%Y%m%d%H%M%S}.nc"
    outp_nc_fn = outp_fn_body+outp_fn_suffix
    outp_nc_fpath = os.path.join(cfg_d["output_dir"], outp_nc_fn)

    output = {}
    output["Global_Attributes"] = ({
                   "full_versionID": product_full_version,
                   "granule_ID": ' ',
                   "provenance": this_pkg_provenance,
                   "file_name": os.path.basename(outp_nc_fpath),
                   "input_product_files": input_product_files,
                   "processing_algorithmID": proc_algID,
                   "UTC_of_file_creation": now_UTC_strrep,
                   "spacecraft_ID": f"PREFIRE0{p_info['sensor_IDval']:1d}",
                   "sensor_ID": p_info["sensor_ID"],
                   "ctime_coverage_start_s": ctime_coverage[0],
                   "ctime_coverage_end_s": ctime_coverage[1],
                   "UTC_coverage_start": UTC_coverage[0],
                   "UTC_coverage_end": UTC_coverage[1],
                   "netCDF_lib_version": netCDF4.getlibversion().split()[0],
                   "orbit_sim_version": p_info["orbit_sim_version"],
                   "SRF_NEdR_version": p_info["SRF_NEdR_version"],
                   "n_usable_calseqs": n_usable_calseqs })

    # Movable-mirror "static" angular offset (with respect to the nadir aperture
    #  boresight):
    # * 1 enc_u = 360./16384. = 0.021972656 degrees (mirror pointing)
    # * since during the slow mirror scan the mirror travels (w.r.t. the
    #    spacecraft body reference frame) from about {a bit +X, +Y} to
    #    {a lot -X, a bit +Y}, a positive encoder unit offset corresponds to a
    #    negative pitch adjustment
    if p_info["sensor_IDval"] == 1:
          # From geostationary radiance comparisons at monitoring sites:
          #   (10./12.)*((0.586 deg/frame)*(10. frames)):
        tmp = 4.88  # [deg]
    else:
          # The -0.4834 deg pitch adjustment (SAT2) is calculated from
          #   * a +22 encoder unit (enc_u) offset -- observed in-orbit -- of the
          # space aperture midpoint, as observed during TIRS2 slow mirror scans
        tmp = -0.4834  # [deg]
          # From geostationary radiance comparisons at monitoring sites:
          #   (10./12.)*((0.586 deg/frame)*(-2.5 frames)):
        tmp += -1.22  # [deg]
    output["TIRS_L0_Pld_Categorization_Group_Attributes"] = ({
                          "mirror_ang_offset0_deg": tmp,
                          "bus_ctime_coverage_start": np.amin(busdat["ctime"]),
                          "bus_ctime_coverage_end": np.amax(busdat["ctime"])})
    output["TIRS_L0_Pld_Categorization"] = o_data

    #== Write to NetCDF-format file:

    mkdir_p(os.path.dirname(outp_nc_fpath))
    
    write_data_fromspec(output, outp_nc_fpath, cfg_d["output_filespecs_fpath"],
                        verbose=True)

    return (n_usable_calseqs, outp_nc_fpath)
