"""
PROC_MODE = 3: Create a Level-1A granule (requires some L0 PROC_MODE 1,2 files)
PROC_MODE = 4: Create part of a Level-1B granule (requires a Level-1A granule)
PROC_MODE = 61: Create a simulated prepped payload tlm file (using the IIM)

This program requires python version 3.6 or later, and is importable as a 
python module.
"""

  # From the Python standard library:
from pathlib import Path
import os
import sys
import argparse
import subprocess

  # From other external Python packages:
import numpy as np

  # Custom utilities:


#--------------------------------------------------------------------------
def main(use_saved_exe, interactive_MatLab):
    """Driver routine."""

    package_top_Path = Path(os.environ["PACKAGE_TOP_DIR"])

    sys.path.append(str(package_top_Path / "source" / "python"))
    from PREFIRE_L1.categorize_L0_payload_tlm import categorize_L0_payload_tlm
    from PREFIRE_L1.curate_L0_input import curate_prepped_L0_input
    from PREFIRE_tools.utils.time import init_leap_s_for_ctimeRefEpoch
    from PREFIRE_tools.utils.unit_conv import hour_to_s

    this_environ = os.environ.copy()

    proc_mode = int(this_environ["PROC_MODE"])

    # Default product_fullver:
    if "PRODUCT_FULLVER" not in this_environ:
        this_environ["PRODUCT_FULLVER"] = "R01_P00"
    elif len(this_environ["PRODUCT_FULLVER"].strip()) == 0:
        this_environ["PRODUCT_FULLVER"] = "R01_P00"

    if proc_mode > 0:   # Requires MatLab
        if proc_mode == 3:
            leap_s_info = init_leap_s_for_ctimeRefEpoch(
                      [2000, 1, 1, 0, 0 ,0], epoch_for_ctime_is_actual_UTC=True)

            pld_tmpcur_fpath, bus_tmpcur_fpath, orb_tmpcur_fpath, ge_t = \
                                            curate_prepped_L0_input(leap_s_info)
            this_environ["PLD_TMPCUR_FPATH"] = pld_tmpcur_fpath
            this_environ["BUS_TMPCUR_FPATH"] = bus_tmpcur_fpath
            this_environ["ORB_TMPCUR_FPATH"] = orb_tmpcur_fpath

               # Set granule_ID from orbit reconstruction data:
            tmp_l = this_environ["GRANULE_START_ID_END"].split('_')
            i_ge = np.argmin(np.abs(np.array(ge_t[0])-float(tmp_l[0])))
            this_environ["GRANULE_START_ID_END"] = '_'.join([tmp_l[0],
                                              f"{ge_t[1][i_ge]:05d}", tmp_l[2]])

            if "PLD_CAT_FPATH" not in this_environ:
                  # Set this cfg dict here since '*_TMPCUR_FPATH' vars cannot be
                  #  conveyed well to `categorize_L0_payload_tlm` otherwise:
                cat_L0p_cfg_d = {}
                cat_L0p_cfg_d["ref_L0_payload_paths"] = (
                                              this_environ["L0_PAYLOAD_FPATHS"])
                cat_L0p_cfg_d["L0_payload_paths"] = (
                                               this_environ["PLD_TMPCUR_FPATH"])
                cat_L0p_cfg_d["ref_L0_bus_paths"] = this_environ["L0_BUS_FPATHS"]
                cat_L0p_cfg_d["L0_bus_paths"] = this_environ["BUS_TMPCUR_FPATH"]
                cat_L0p_cfg_d["ref_L0_orbit_paths"] = (
                                                this_environ["L0_ORBIT_FPATHS"])
                cat_L0p_cfg_d["L0_orbit_paths"] = (
                                               this_environ["ORB_TMPCUR_FPATH"])
                cat_L0p_cfg_d["output_filespecs_fpath"] = os.path.join(
                                             this_environ["ANCILLARY_DATA_DIR"],
                                         "categorization_L0_pld_filespecs.json")
                cat_L0p_cfg_d["output_dir"] = this_environ["OUTPUT_DIR"]
                cat_L0p_cfg_d["product_fullver"] = (
                                                this_environ["PRODUCT_FULLVER"])

                tmp_l = this_environ["GRANULE_START_ID_END"].split('_')
                cat_L0p_cfg_d["ctime_range_to_proc"] = (
                                                float(tmp_l[0])-16.99*hour_to_s,
                                                float(tmp_l[2])+1.59*hour_to_s)

                n_calseqs, pld_cat_fpath = categorize_L0_payload_tlm(leap_s_info,
                                                                  cat_L0p_cfg_d)
                if abs(n_calseqs) < 3:
                    msg = ("Too few ({}; minimum is 3) calibration "
                           "sequences were found.".format(n_calseqs))
                    if use_saved_exe:  # Operational case
                        raise RuntimeError(msg)
                    else:
                        this_environ["THERE_IS_A_PRE_ERROR"] = "yes"
                        print(msg)  # do not shut down the pipeline
                this_environ["PLD_CAT_FPATH"] = pld_cat_fpath

        OS_is_MSWin = ( not os.sep == '/' )
          # Determine separator for use in a search PATH:
        if OS_is_MSWin:
            Psep = ';'
        else:
            Psep = ':'

        MatLab_src_dir = str(package_top_Path / "source" / "matlab" /
                             "functions")
        tools_src_dir = str(package_top_Path / "source" / "matlab" /
                            "PREFIRE_tools" / "functions")
        tools_anc_dir0 = str(package_top_Path / "source" / "matlab" /
                             "PREFIRE_tools" / ".." / ".." / "dist" /
                             "ancillary")
        this_environ["TOOLS_ANC_DIR"] = os.path.realpath(tools_anc_dir0)
        app_MatLab_prefix = "app_matlab"

        dist_Path = package_top_Path / "dist"

        if interactive_MatLab:
            this_environ["MATLABPATH"] = (
                     f"{MatLab_src_dir}{Psep}{dist_Path}{Psep}{tools_src_dir}")
            cmd = ["matlab", "-singleCompThread"]
        else:
            if use_saved_exe:
                cmd = [str(dist_Path / f"{app_MatLab_prefix}_run")]
            else:
                this_environ["MATLABPATH"] = (
                     f"{MatLab_src_dir}{Psep}{dist_Path}{Psep}{tools_src_dir}")
                if OS_is_MSWin:
                    cmd = ["matlab", "-singleCompThread", "-noFigureWindows", "-batch",
                           f"{app_MatLab_prefix}_run"]
                else:
                    cmd = ["matlab", "-singleCompThread", "-nodisplay", "-batch",
                           f"{app_MatLab_prefix}_run"]

        try:
            result = subprocess.run(cmd, env=this_environ, check=True)
        except subprocess.CalledProcessError as e:
            if use_saved_exe:  # Operational case
                raise e from None
            else:
                pass  # Error info should have been printed by MatLab app,
                      #  do not shut down the pipeline


if __name__ == "__main__":
    # Process arguments:
    arg_description = ("PROC_MODE = 3: Create a Level-1A granule\n"
                       "PROC_MODE = 4: Create all/part of a Level-1B granule")
    arg_parser = argparse.ArgumentParser(description=arg_description)
    arg_parser.add_argument("-s", "--use_saved_exe", action="store_true",
                            help="Use the MatLab Runtime to execute routines "
                                 "stored in a Matlab executable file.")
    arg_parser.add_argument("-i", "--interactive-no_display",
                            dest="interactive", action="store_true",
                            help="Set environment, then run MatLab in "
                               "interactive mode (but with no fancy display).") 

    args = arg_parser.parse_args()

    # Check for argument sanity:
    if args.use_saved_exe and args.interactive:
        raise ValueError("Arguments -s and -i cannot be used together.")

    # Run driver:
    main(args.use_saved_exe, args.interactive)
