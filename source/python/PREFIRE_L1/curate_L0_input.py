"""
From one or more prepped Level-0 file(s), create a curated version (i.e., all
the frames chronologically ordered and de-duplicated), and write the resulting
information to a new NetCDF-format file.  Do this for both payload and bus
telemetry information (and orbit reconstruction data).

This program requires Python version 3.6 or later, and is importable as a
Python module.
"""

  # From the Python standard library:
from pathlib import Path
import os

  # From other external Python packages:
import netCDF4

  # Custom utilities:
from PREFIRE_PRD_GEN.file_concatenation import dedup_and_concat_along_tdim


#--------------------------------------------------------------------------
def curate_prepped_L0_input(leap_s_info, ext_cfg_d=None):
    """Driver routine."""

    cfg_d = {}
    if ext_cfg_d is None:
        cfg_d["L0_payload_paths"] = os.environ["L0_PAYLOAD_FPATHS"]
        cfg_d["L0_bus_paths"] = os.environ["L0_BUS_FPATHS"]
        cfg_d["L0_orbit_paths"] = os.environ["L0_ORBIT_FPATHS"]
        cfg_d["output_dir"] = os.environ["OUTPUT_DIR"]
    else:
        cfg_d = ext_cfg_d

    outp_L0pld_fpath = None
    outp_L0bus_fpath = None

    # Curate any prepped payload telemetry file input:
    if len(cfg_d["L0_payload_paths"].strip()) > 0:
        input_p_fpaths = [x.strip() for x in
                          cfg_d["L0_payload_paths"].split("||")]
        outp_L0pld_fpath = dedup_and_concat_along_tdim(input_p_fpaths,
                                     Path(cfg_d["output_dir"]), "atrack_scipkt",
                                   ("ctime", "TIRS_L0_Data", None), leap_s_info)

    # Curate any prepped bus telemetry file input:
    if len(cfg_d["L0_bus_paths"].strip()) > 0:
        input_b_fpaths = [x.strip() for x in cfg_d["L0_bus_paths"].split("||")]
        outp_L0bus_fpath = dedup_and_concat_along_tdim(input_b_fpaths,
                                           Path(cfg_d["output_dir"]), "atrack",
                                  ("ctime", "SCbus_L0_Data", None), leap_s_info)
        ge_cfg = (input_b_fpaths, "SCbus_L0_Data")

    # Curate any prepped orbit reconstruction file input:
    if len(cfg_d["L0_orbit_paths"].strip()) > 0:
        input_b_fpaths = [x.strip() for x in cfg_d["L0_orbit_paths"].split("||")]
        outp_L0orb_fpath = dedup_and_concat_along_tdim(input_b_fpaths,
                                           Path(cfg_d["output_dir"]), "atrack",
                                  ("ctime", "Orbit_L0_Data", None), leap_s_info)
        ge_cfg = (input_b_fpaths, "Orbit_L0_Data")

    ge_ct_l = []
    ge_ID_l = []
    for fpath in ge_cfg[0]:
            with netCDF4.Dataset(fpath, 'r') as ds:
                ge_ct_l.extend(ds.groups[ge_cfg[1]].
                               variables["granule_edge_ctime"][...])
                ge_ID_l.extend(ds.groups[ge_cfg[1]].
                               variables["granule_ID"][...])

    return (outp_L0pld_fpath, outp_L0bus_fpath, outp_L0orb_fpath,
            (ge_ct_l, ge_ID_l))
