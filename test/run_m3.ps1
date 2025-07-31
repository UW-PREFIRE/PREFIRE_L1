## IMPORTANT: Only run this script from the directory it resides in, i.e. with
##             .\run_m3.ps1 [-i]

##===========================================================================##
## This script contains hardwired information necessary for this algorithm's
##  delivery to and testing within the SDPS (Science Data Processing System).
##
## ** In general, do not push changes to this file to its primary git
##     repository (exceptions include adding a new environment var for
##     algorithm config) **
##
## ++ Instead, make a LOCAL copy of this script (e.g., my_run_m3.ps1; do not
##     push that local copy to the primary git repository either) and modify
##     and run that for general algorithm testing and development.
##===========================================================================##

param(
  [Parameter()]
  [switch]$i
)

# Determine the absolute path of the current working directory:
#  (this is typically the package test/ directory)
$base_dir = $pwd

# NOTE: Set the input/output directories to absolute paths (relative to the
#        current working directory, 'base_dir').

$L0_dir = "$base_dir\inputs"

$L0_granule_cfg_str1 = "782505241.569744_00000_782510951.185644:$L0_dir\prefire_01_payload_tlm_20241017000000_20241017235959_20241210183433.nc,$L0_dir\prefire_01_bus_tlm_20241017010000_20241017015959_20241017041244.nc||$L0_dir\prefire_01_bus_tlm_20241017020000_20241017025959_20241017081059.nc||$L0_dir\prefire_01_bus_tlm_20241017030000_20241017035959_20241017081059.nc||$L0_dir\prefire_01_bus_tlm_20241017040000_20241017045959_20241017161401.nc||$L0_dir\prefire_01_bus_tlm_20241017050000_20241017055959_20241017161401.nc||$L0_dir\prefire_01_bus_tlm_20241017060000_20241017065959_20241017161401.nc||$L0_dir\prefire_01_bus_tlm_20241017070000_20241017075959_20241017161422.nc||$L0_dir\prefire_01_bus_tlm_20241017080000_20241017085959_20241017161422.nc||$L0_dir\prefire_01_bus_tlm_20241017090000_20241017095959_20241017161422.nc||$L0_dir\prefire_01_bus_tlm_20241017100000_20241017105959_20241017161422.nc||$L0_dir\prefire_01_bus_tlm_20241017110000_20241017115959_20241017161422.nc||$L0_dir\prefire_01_bus_tlm_20241017120000_20241017125959_20241017161422.nc||$L0_dir\prefire_01_bus_tlm_20241017130000_20241017135959_20241017201328.nc||$L0_dir\prefire_01_bus_tlm_20241017140000_20241017145959_20241017201328.nc||$L0_dir\prefire_01_bus_tlm_20241017150000_20241017155959_20241017201328.nc||$L0_dir\prefire_01_bus_tlm_20241017160000_20241017165959_20241017201328.nc||$L0_dir\prefire_01_bus_tlm_20241017170000_20241017175959_20241018121014.nc||$L0_dir\prefire_01_bus_tlm_20241017180000_20241017185959_20241018121014.nc||$L0_dir\prefire_01_bus_tlm_20241017190000_20241017195959_20241018121014.nc||$L0_dir\prefire_01_bus_tlm_20241017200000_20241017205959_20241018121014.nc||$L0_dir\prefire_01_bus_tlm_20241017210000_20241017215959_20241018121014.nc,$L0_dir\prefire_01_orbit_reconst_20241017010000_20241017015959_20250116193755.nc||$L0_dir\prefire_01_orbit_reconst_20241017020000_20241017025959_20250116193755.nc||$L0_dir\prefire_01_orbit_reconst_20241017030000_20241017035959_20241017081059.nc||$L0_dir\prefire_01_orbit_reconst_20241017040000_20241017045959_20241017161401.nc||$L0_dir\prefire_01_orbit_reconst_20241017050000_20241017055959_20241017161401.nc||$L0_dir\prefire_01_orbit_reconst_20241017060000_20241017065959_20241017161401.nc||$L0_dir\prefire_01_orbit_reconst_20241017070000_20241017075959_20241017161422.nc||$L0_dir\prefire_01_orbit_reconst_20241017080000_20241017085959_20241017161422.nc||$L0_dir\prefire_01_orbit_reconst_20241017090000_20241017095959_20241017161422.nc||$L0_dir\prefire_01_orbit_reconst_20241017100000_20241017105959_20241017161422.nc||$L0_dir\prefire_01_orbit_reconst_20241017110000_20241017115959_20241017161422.nc||$L0_dir\prefire_01_orbit_reconst_20241017120000_20241017125959_20241017161422.nc||$L0_dir\prefire_01_orbit_reconst_20241017130000_20241017135959_20241017201328.nc||$L0_dir\prefire_01_orbit_reconst_20241017140000_20241017145959_20241017201328.nc||$L0_dir\prefire_01_orbit_reconst_20241017150000_20241017155959_20241017201328.nc||$L0_dir\prefire_01_orbit_reconst_20241017160000_20241017165959_20241017201328.nc||$L0_dir\prefire_01_orbit_reconst_20241017170000_20241017175959_20241018121014.nc||$L0_dir\prefire_01_orbit_reconst_20241017180000_20241017185959_20241018121014.nc||$L0_dir\prefire_01_orbit_reconst_20241017190000_20241017195959_20241018121014.nc||$L0_dir\prefire_01_orbit_reconst_20241017200000_20241017205959_20241018121014.nc||$L0_dir\prefire_01_orbit_reconst_20241017210000_20241017215959_20241018121014.nc"

$L0_granule_cfg_str2 = "781811871.951264_00000_781817583.267378:$L0_dir\prefire_02_payload_tlm_20241009000000_20241009235959_20241210183452.nc,$L0_dir\prefire_02_bus_tlm_20241009000000_20241009005959_20241009080503.nc||$L0_dir\prefire_02_payload_tlm_20241009000000_20241009235959_20241210183452.nc,$L0_dir\prefire_02_bus_tlm_20241009020000_20241009025959_20241009080503.nc||$L0_dir\prefire_02_bus_tlm_20241009030000_20241009035959_20241009080503.nc||$L0_dir\prefire_02_bus_tlm_20241009040000_20241009045959_20241009080503.nc||$L0_dir\prefire_02_bus_tlm_20241009050000_20241009055959_20241009080503.nc||$L0_dir\prefire_02_bus_tlm_20241009060000_20241009065959_20241009120446.nc||$L0_dir\prefire_02_bus_tlm_20241009070000_20241009075959_20241009120446.nc||$L0_dir\prefire_02_bus_tlm_20241009080000_20241009085959_20241010000530.nc||$L0_dir\prefire_02_bus_tlm_20241009090000_20241009095959_20241010000530.nc||$L0_dir\prefire_02_bus_tlm_20241009100000_20241009105959_20241010000530.nc||$L0_dir\prefire_02_bus_tlm_20241009110000_20241009115959_20241010000530.nc||$L0_dir\prefire_02_bus_tlm_20241009120000_20241009125959_20241010000530.nc||$L0_dir\prefire_02_bus_tlm_20241009130000_20241009135959_20241010000530.nc||$L0_dir\prefire_02_bus_tlm_20241009140000_20241009145959_20241010000530.nc||$L0_dir\prefire_02_bus_tlm_20241009150000_20241009155959_20241010000544.nc||$L0_dir\prefire_02_bus_tlm_20241009160000_20241009165959_20241010000544.nc||$L0_dir\prefire_02_bus_tlm_20241009170000_20241009175959_20241010000544.nc||$L0_dir\prefire_02_bus_tlm_20241009180000_20241009185959_20241010120534.nc||$L0_dir\prefire_02_bus_tlm_20241009190000_20241009195959_20241010120534.nc||$L0_dir\prefire_02_bus_tlm_20241009200000_20241009205959_20241010120534.nc||$L0_dir\prefire_02_bus_tlm_20241009210000_20241009215959_20241010120534.nc,$L0_dir\prefire_02_orbit_reconst_20241009000000_20241009005959_20250108202214.nc||$L0_dir\prefire_02_orbit_reconst_20241009010000_20241009015959_20250108202214.nc||$L0_dir\prefire_02_orbit_reconst_20241009020000_20241009025959_20241009080503.nc||$L0_dir\prefire_02_orbit_reconst_20241009030000_20241009035959_20241009080503.nc||$L0_dir\prefire_02_orbit_reconst_20241009040000_20241009045959_20241009080503.nc||$L0_dir\prefire_02_orbit_reconst_20241009050000_20241009055959_20241009080503.nc||$L0_dir\prefire_02_orbit_reconst_20241009060000_20241009065959_20241009120446.nc||$L0_dir\prefire_02_orbit_reconst_20241009070000_20241009075959_20241009120446.nc||$L0_dir\prefire_02_orbit_reconst_20241009080000_20241009085959_20241010000530.nc||$L0_dir\prefire_02_orbit_reconst_20241009090000_20241009095959_20241010000530.nc||$L0_dir\prefire_02_orbit_reconst_20241009100000_20241009105959_20241010000530.nc||$L0_dir\prefire_02_orbit_reconst_20241009110000_20241009115959_20241010000530.nc||$L0_dir\prefire_02_orbit_reconst_20241009120000_20241009125959_20241010000530.nc||$L0_dir\prefire_02_orbit_reconst_20241009130000_20241009135959_20241010000530.nc||$L0_dir\prefire_02_orbit_reconst_20241009140000_20241009145959_20241010000530.nc||$L0_dir\prefire_02_orbit_reconst_20241009150000_20241009155959_20241010000544.nc||$L0_dir\prefire_02_orbit_reconst_20241009160000_20241009165959_20241010000544.nc||$L0_dir\prefire_02_orbit_reconst_20241009170000_20241009175959_20241010000544.nc||$L0_dir\prefire_02_orbit_reconst_20241009180000_20241009185959_20241010120534.nc||$L0_dir\prefire_02_orbit_reconst_20241009190000_20241009195959_20241010120534.nc||$L0_dir\prefire_02_orbit_reconst_20241009200000_20241009205959_20241010120534.nc||$L0_dir\prefire_02_orbit_reconst_20241009210000_20241009215959_20241010120534.nc"


# Specify that numpy, scipy, et cetera should not use more than one thread or
#  process):
$env:MKL_NUM_THREADS = '1'
$env:NUMEXPR_NUM_THREADS = '1'
$env:OMP_NUM_THREADS = '1'
$env:VECLIB_MAXIMUM_THREADS = '1'
$env:OPENBLAS_NUM_THREADS = '1'

# Some environment vars that convey configuration info to the algorithm:

$this_top_dir = [IO.Path]::GetFullPath("$base_dir\..")

$env:PACKAGE_TOP_DIR = "$this_top_dir"
$env:ANCILLARY_DATA_DIR = "$this_top_dir\dist\ancillary"
$env:INSTRUMENT_MODEL_DIR = "$env:ANCILLARY_DATA_DIR\instrument_model"

#====== For certain types of development and testing ======
#$env:SRF_DISAMBIG_STR = "SRF_v12_2024-04-18"
#==========================================================

$env:OUTPUT_DIR = "$base_dir\outputs\m3"

$env:PROC_MODE = '3'

  # * Only increment 'Rxx' when the resulting products will be DAAC-ingested
$env:PRODUCT_FULLVER = "R01_P00"
  # Special form ('R00_Syy') when processing simulated observations:
#$env:PRODUCT_FULLVER = "R00_S01"

#= Processing mode #3: Create a Level-1A granule (requires some PROC_MODE 1,2
#                      files).
#  NOTE: Glue code must determine and set the values of GRANULE_START_ID_END,
#        L0_PAYLOAD_FPATHS, and L0_BUS_FPATHS and make sure the *_FPATHS files
#        are copied to the working directory.
#
# for GRANULE_START_ID_END:
#   composed of 3 numeric strings separated by underscores:
#    1) ctime of granule start (edge) [microseconds]
#    2) granule ID, 5 chars in length, left-padded with zeros
#    3) ctime of granule end (edge) [microseconds]
#
# for *_FPATHS:
#   can contain more than one filepath, in which case they should be
#    separated with '||'

# Check if output file directory exists; if not, bail:
$tmpdir = "$env:OUTPUT_DIR"
If (-not (Test-Path -Path $tmpdir)) {
  throw "Output directory does not exist: $tmpdir"
}

# Execute script that writes a new 'prdgit_version.txt', which contains
#  product moniker(s) and current (latest) git hash(es) that are part of the
#  provenance of this package's product(s).
# *** This step should not be done within the SDPS, since that file is
#     created just before delivery to the SDPS.
If (-not (Test-Path -Path "$this_top_dir\dist\for_SDPS_delivery.txt")) {
  python "$this_top_dir\dist\determine_prdgit.py"
}

$items = @($L0_granule_cfg_str1, $L0_granule_cfg_str2)
foreach ($cfg_str in $items) {
  $split_strarr = $cfg_str -split ";"

  $env:GRANULE_START_ID_END = $split_strarr[0]

  $tmp_strarr = $split_strarr[1] -split ","
  $env:L0_PAYLOAD_FPATHS = $tmp_strarr[0]
  $env:L0_BUS_FPATHS = $tmp_strarr[1]
  $env:L0_ORBIT_FPATHS = $tmp_strarr[2]

  if ($i.IsPresent) {
    python "$this_top_dir\dist\produce_L1.py" -i
  } else {
    python "$this_top_dir\dist\produce_L1.py"
  }
}
