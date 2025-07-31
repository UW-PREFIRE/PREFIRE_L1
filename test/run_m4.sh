#!/usr/bin/env bash

## IMPORTANT: Only run this script from the directory it resides in, i.e. with
##             ./run_m4.sh    OR    bash run_m4.sh

##===========================================================================##
## This script contains hardwired information necessary for this algorithm's
##  delivery to and testing within the SDPS (Science Data Processing System).
##
## ** In general, do not push changes to this file to its primary git
##     repository (exceptions include adding a new environment var for
##     algorithm config) **
##
## ++ Instead, make a LOCAL copy of this script (e.g., my_run_m4.sh; do not
##     push that local copy to the primary git repository either) and modify
##     and run that for general algorithm testing and development.
##===========================================================================##

absfpath() {
  # Generate absolute filepath from a relative (or even an absolute) filepath.
  #
  # Based on (circa Oct 2023) https://stackoverflow.com/questions/3915040/how-to-obtain-the-absolute-path-of-a-file-via-shell-bash-zsh-sh
  # 
  # $1     : a relative (or even an absolute) filepath
  # Returns the corresponding absolute filepath.
  if [ -d "$1" ]; then
    # dir
    (cd "$1"; pwd)
  elif [ -f "$1" ]; then
    # file
    if [[ $1 = /* ]]; then
      echo "$1"
    elif [[ $1 == */* ]]; then
      echo "$(cd "${1%/*}"; pwd)/${1##*/}"
    else
      echo "$(pwd)/$1"
    fi
  fi
}

activate_conda_env () {
  . "$1"/bin/activate;
}

deactivate_conda_env () {
  . "$1"/bin/deactivate;
}

#set -ve;  # Exit on the first error, and print out commands as we execute them
set -e;  # Exit on the first error

# Determine the absolute path of the current working directory:
#  (this is typically the package test/ directory)
readonly base_dir="$(absfpath ".")";

hn=`hostname -s`;  # Hostname

# NOTE: Set the input/output directories to absolute paths (relative to the
#        current working directory, 'base_dir').

non_SDPS_hostname="longwave";

#L0_dir="${base_dir}/inputs";
L0_dir="${base_dir}/outputs/m3";

L1A_granule_cfg_str1="ATRACK_IDXRANGE_0BASED_INCLUSIVE:4000:4100,${L0_dir}/raw-PREFIRE_SAT1_1A-RAD_R01_P00_20241017183356_02035.nc";
L1A_granule_cfg_str2="ATRACK_IDXRANGE_0BASED_INCLUSIVE:400:500,${L0_dir}/raw-PREFIRE_SAT2_1A-RAD_R01_P00_20241009175747_02077.nc";


# Specify that numpy, scipy, et cetera should not use more than one thread or
#  process):
MKL_NUM_THREADS=1;
NUMEXPR_NUM_THREADS=1;
OMP_NUM_THREADS=1;
VECLIB_MAXIMUM_THREADS=1;
OPENBLAS_NUM_THREADS=1;
export MKL_NUM_THREADS NUMEXPR_NUM_THREADS OMP_NUM_THREADS;
export VECLIB_MAXIMUM_THREADS OPENBLAS_NUM_THREADS;

# Some environment vars that convey configuration info to the algorithm:

this_top_dir="$(absfpath "${base_dir}/..")";

PACKAGE_TOP_DIR="${this_top_dir}";
ANCILLARY_DATA_DIR="${this_top_dir}/dist/ancillary";
INSTRUMENT_MODEL_DIR="${ANCILLARY_DATA_DIR}/instrument_model";

DEM_ROOT_DIR=/data/users/mmm/DEM/copernicus-dem-90m/tiles;

OUTPUT_DIR="${base_dir}/outputs/m4";

PROC_MODE=4;

  # * Only increment 'Rxx' when the resulting products will be DAAC-ingested
PRODUCT_FULLVER="R01_P00";
  # Special form ('R00_Syy') when processing simulated observations:
#PRODUCT_FULLVER="R00_S01";

export PACKAGE_TOP_DIR ANCILLARY_DATA_DIR INSTRUMENT_MODEL_DIR DEM_ROOT_DIR;
export OUTPUT_DIR PROC_MODE PRODUCT_FULLVER;

#= Processing mode #4: Create part of a Level-1B granule (requires a Level-1A
#                      granule).
#
# for the ATRACK_IDX_RANGE_0BI string:
#   The last two colon-separated items describe the frame subset of this
#    granule to process -- they are the 0-based inclusive frame indices of the
#    desired frame subset, where the end index can also be set to 'END'.

# Check if output file directory exists; if not, bail:
tmpdir="${OUTPUT_DIR}";
test -d "${tmpdir}" || { echo "Output directory does not exist: ${tmpdir}"; exit 1; }

# If custom conda environment files exist, activate that conda environment:
conda_env_dir="${this_top_dir}/dist/c_env_for_PREFIRE_L1";
if [ -d "${conda_env_dir}" ]; then
   activate_conda_env "${conda_env_dir}";
fi

# Execute script that writes a new 'prdgit_version.txt', which contains
#  product moniker(s) and current (latest) git hash(es) that are part of the
#  provenance of this package's product(s).
# *** This step should not be done within the SDPS, since that file is
#     created just before delivery to the SDPS.
if [ ! -f "${this_top_dir}/dist/for_SDPS_delivery.txt" ]; then
   python "${this_top_dir}/dist/determine_prdgit.py";
fi

# Execute any necessary machine setup instructions ((un)loading modules, etc.):
if [ "x$hn" = "x$non_SDPS_hostname" ]; then
   . "${this_top_dir}/dist/perform_machine_setup.sh";
fi

for cfg_str in ${L1A_granule_cfg_str1} ${L1A_granule_cfg_str2}
do
   ATRACK_IDX_RANGE_0BI=${cfg_str%,*};
   L1A_RAD_FILE=${cfg_str##*,};

   export L1A_RAD_FILE ATRACK_IDX_RANGE_0BI;

   if [ "x$1" = "x-i" ]; then
      python "${this_top_dir}/dist/produce_L1.py" -i;
   else
      if [ "x$hn" = "x$non_SDPS_hostname" ]; then
         #!!! The MATLAB module on longwave does not seem to set these needed
         #    paths (for standalone M-executable execution) in any way, so do
         #    it here:
         MATLAB_BASEDIR=/opt/matlab/r2021b;
         LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$MATLAB_BASEDIR/bin/glnxa64:$MATLAB_BASEDIR/runtime/glnxa64";
         export LD_LIBRARY_PATH;
      fi

      python "${this_top_dir}/dist/produce_L1.py" $@;
   fi
done

# If custom conda environment files exist, DEactivate that conda environment:
if [ -d "${conda_env_dir}" ]; then
   deactivate_conda_env "${conda_env_dir}";
fi
