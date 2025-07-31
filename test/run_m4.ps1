## IMPORTANT: Only run this script from the directory it resides in, i.e. with
##             .\run_m4.ps1

##===========================================================================##
## This script contains hardwired information necessary for this algorithm's
##  delivery to and testing within the SDPS (Science Data Processing System).
##
## ** In general, do not push changes to this file to its primary git
##     repository (exceptions include adding a new environment var for
##     algorithm config) **
##
## ++ Instead, make a LOCAL copy of this script (e.g., my_run_m4.ps1; do not
##     push that local copy to the primary git repository either) and modify
##     and run that for general algorithm testing and development.
##===========================================================================##

# Determine the absolute path of the current working directory:
#  (this is typically the package test/ directory)
$base_dir = $pwd

# NOTE: Set the input/output directories to absolute paths (relative to the
#        current working directory, 'base_dir').

$L0_dir = "$base_dir\inputs"
L0_dir = "$base_dir\outputs\m3"

$L1A_granule_cfg_str1 = "ATRACK_IDXRANGE_0BASED_INCLUSIVE:4000:4100,$L0_dir\raw-PREFIRE_SAT1_1A-RAD_R01_P00_20241017183356_02035.nc"
$L1A_granule_cfg_str2 = "ATRACK_IDXRANGE_0BASED_INCLUSIVE:400:500,$L0_dir\raw-PREFIRE_SAT2_1A-RAD_R01_P00_20241009175747_02077.nc"


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

$env:DEM_ROOT_DIR = "C:\Users\mmm\DEM\copernicus-dem-90m\tiles"

$env:OUTPUT_DIR = "$base_dir\outputs\m4"

$env:PROC_MODE = '4'

  # * Only increment 'Rxx' when the resulting products will be DAAC-ingested
$env:PRODUCT_FULLVER = "R01_P00"
  # Special form ('R00_Syy') when processing simulated observations:
#$env:PRODUCT_FULLVER = "R00_S01"

#= Processing mode #4: Create part of a Level-1B granule (requires a Level-1A
#                      granule).
#
# for the ATRACK_IDX_RANGE_0BI string:
#   The last two colon-separated items describe the frame subset of this
#    granule to process -- they are the 0-based inclusive frame indices of the
#    desired frame subset, where the end index can also be set to 'END'.

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

$items = @($L1A_granule_cfg_str1, $L1A_granule_cfg_str2)
foreach ($cfg_str in $items) {
  $split_strarr = $cfg_str -split ","

  $env:ATRACK_IDX_RANGE_0BI = $split_strarr[0]
  $env:L1A_RAD_FILE = $split_strarr[1]

  #python "$this_top_dir\dist\produce_L1.py" -i
  python "$this_top_dir\dist\produce_L1.py"
}
