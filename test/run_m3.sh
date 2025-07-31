#!/usr/bin/env bash

## IMPORTANT: Only run this script from the directory it resides in, i.e. with
##             ./run_m3.sh    OR    bash run_m3.sh

##===========================================================================##
## This script contains hardwired information necessary for this algorithm's
##  delivery to and testing within the SDPS (Science Data Processing System).
##
## ** In general, do not push changes to this file to its primary git
##     repository (exceptions include adding a new environment var for
##     algorithm config) **
##
## ++ Instead, make a LOCAL copy of this script (e.g., my_run_m3.sh; do not
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

L0_dir="${base_dir}/inputs";

L0_granule_cfg_str1="782505241.569744_00000_782510951.185644:${L0_dir}/prefire_01_payload_tlm_20241017000000_20241017235959_20241210183433.nc,${L0_dir}/prefire_01_bus_tlm_20241017010000_20241017015959_20241017041244.nc||${L0_dir}/prefire_01_bus_tlm_20241017020000_20241017025959_20241017081059.nc||${L0_dir}/prefire_01_bus_tlm_20241017030000_20241017035959_20241017081059.nc||${L0_dir}/prefire_01_bus_tlm_20241017040000_20241017045959_20241017161401.nc||${L0_dir}/prefire_01_bus_tlm_20241017050000_20241017055959_20241017161401.nc||${L0_dir}/prefire_01_bus_tlm_20241017060000_20241017065959_20241017161401.nc||${L0_dir}/prefire_01_bus_tlm_20241017070000_20241017075959_20241017161422.nc||${L0_dir}/prefire_01_bus_tlm_20241017080000_20241017085959_20241017161422.nc||${L0_dir}/prefire_01_bus_tlm_20241017090000_20241017095959_20241017161422.nc||${L0_dir}/prefire_01_bus_tlm_20241017100000_20241017105959_20241017161422.nc||${L0_dir}/prefire_01_bus_tlm_20241017110000_20241017115959_20241017161422.nc||${L0_dir}/prefire_01_bus_tlm_20241017120000_20241017125959_20241017161422.nc||${L0_dir}/prefire_01_bus_tlm_20241017130000_20241017135959_20241017201328.nc||${L0_dir}/prefire_01_bus_tlm_20241017140000_20241017145959_20241017201328.nc||${L0_dir}/prefire_01_bus_tlm_20241017150000_20241017155959_20241017201328.nc||${L0_dir}/prefire_01_bus_tlm_20241017160000_20241017165959_20241017201328.nc||${L0_dir}/prefire_01_bus_tlm_20241017170000_20241017175959_20241018121014.nc||${L0_dir}/prefire_01_bus_tlm_20241017180000_20241017185959_20241018121014.nc||${L0_dir}/prefire_01_bus_tlm_20241017190000_20241017195959_20241018121014.nc||${L0_dir}/prefire_01_bus_tlm_20241017200000_20241017205959_20241018121014.nc||${L0_dir}/prefire_01_bus_tlm_20241017210000_20241017215959_20241018121014.nc,${L0_dir}/prefire_01_orbit_reconst_20241017010000_20241017015959_20250116193755.nc||${L0_dir}/prefire_01_orbit_reconst_20241017020000_20241017025959_20250116193755.nc||${L0_dir}/prefire_01_orbit_reconst_20241017030000_20241017035959_20241017081059.nc||${L0_dir}/prefire_01_orbit_reconst_20241017040000_20241017045959_20241017161401.nc||${L0_dir}/prefire_01_orbit_reconst_20241017050000_20241017055959_20241017161401.nc||${L0_dir}/prefire_01_orbit_reconst_20241017060000_20241017065959_20241017161401.nc||${L0_dir}/prefire_01_orbit_reconst_20241017070000_20241017075959_20241017161422.nc||${L0_dir}/prefire_01_orbit_reconst_20241017080000_20241017085959_20241017161422.nc||${L0_dir}/prefire_01_orbit_reconst_20241017090000_20241017095959_20241017161422.nc||${L0_dir}/prefire_01_orbit_reconst_20241017100000_20241017105959_20241017161422.nc||${L0_dir}/prefire_01_orbit_reconst_20241017110000_20241017115959_20241017161422.nc||${L0_dir}/prefire_01_orbit_reconst_20241017120000_20241017125959_20241017161422.nc||${L0_dir}/prefire_01_orbit_reconst_20241017130000_20241017135959_20241017201328.nc||${L0_dir}/prefire_01_orbit_reconst_20241017140000_20241017145959_20241017201328.nc||${L0_dir}/prefire_01_orbit_reconst_20241017150000_20241017155959_20241017201328.nc||${L0_dir}/prefire_01_orbit_reconst_20241017160000_20241017165959_20241017201328.nc||${L0_dir}/prefire_01_orbit_reconst_20241017170000_20241017175959_20241018121014.nc||${L0_dir}/prefire_01_orbit_reconst_20241017180000_20241017185959_20241018121014.nc||${L0_dir}/prefire_01_orbit_reconst_20241017190000_20241017195959_20241018121014.nc||${L0_dir}/prefire_01_orbit_reconst_20241017200000_20241017205959_20241018121014.nc||${L0_dir}/prefire_01_orbit_reconst_20241017210000_20241017215959_20241018121014.nc";

L0_granule_cfg_str2="781811871.951264_00000_781817583.267378:${L0_dir}/prefire_02_payload_tlm_20241009000000_20241009235959_20241210183452.nc,${L0_dir}/prefire_02_bus_tlm_20241009000000_20241009005959_20241009080503.nc||${L0_dir}/prefire_02_bus_tlm_20241009010000_20241009015959_20241009080503.nc||${L0_dir}/prefire_02_bus_tlm_20241009020000_20241009025959_20241009080503.nc||${L0_dir}/prefire_02_bus_tlm_20241009030000_20241009035959_20241009080503.nc||${L0_dir}/prefire_02_bus_tlm_20241009040000_20241009045959_20241009080503.nc||${L0_dir}/prefire_02_bus_tlm_20241009050000_20241009055959_20241009080503.nc||${L0_dir}/prefire_02_bus_tlm_20241009060000_20241009065959_20241009120446.nc||${L0_dir}/prefire_02_bus_tlm_20241009070000_20241009075959_20241009120446.nc||${L0_dir}/prefire_02_bus_tlm_20241009080000_20241009085959_20241010000530.nc||${L0_dir}/prefire_02_bus_tlm_20241009090000_20241009095959_20241010000530.nc||${L0_dir}/prefire_02_bus_tlm_20241009100000_20241009105959_20241010000530.nc||${L0_dir}/prefire_02_bus_tlm_20241009110000_20241009115959_20241010000530.nc||${L0_dir}/prefire_02_bus_tlm_20241009120000_20241009125959_20241010000530.nc||${L0_dir}/prefire_02_bus_tlm_20241009130000_20241009135959_20241010000530.nc||${L0_dir}/prefire_02_bus_tlm_20241009140000_20241009145959_20241010000530.nc||${L0_dir}/prefire_02_bus_tlm_20241009150000_20241009155959_20241010000544.nc||${L0_dir}/prefire_02_bus_tlm_20241009160000_20241009165959_20241010000544.nc||${L0_dir}/prefire_02_bus_tlm_20241009170000_20241009175959_20241010000544.nc||${L0_dir}/prefire_02_bus_tlm_20241009180000_20241009185959_20241010120534.nc||${L0_dir}/prefire_02_bus_tlm_20241009190000_20241009195959_20241010120534.nc||${L0_dir}/prefire_02_bus_tlm_20241009200000_20241009205959_20241010120534.nc||${L0_dir}/prefire_02_bus_tlm_20241009210000_20241009215959_20241010120534.nc,${L0_dir}/prefire_02_orbit_reconst_20241009000000_20241009005959_20250108202214.nc||${L0_dir}/prefire_02_orbit_reconst_20241009010000_20241009015959_20250108202214.nc||${L0_dir}/prefire_02_orbit_reconst_20241009020000_20241009025959_20241009080503.nc||${L0_dir}/prefire_02_orbit_reconst_20241009030000_20241009035959_20241009080503.nc||${L0_dir}/prefire_02_orbit_reconst_20241009040000_20241009045959_20241009080503.nc||${L0_dir}/prefire_02_orbit_reconst_20241009050000_20241009055959_20241009080503.nc||${L0_dir}/prefire_02_orbit_reconst_20241009060000_20241009065959_20241009120446.nc||${L0_dir}/prefire_02_orbit_reconst_20241009070000_20241009075959_20241009120446.nc||${L0_dir}/prefire_02_orbit_reconst_20241009080000_20241009085959_20241010000530.nc||${L0_dir}/prefire_02_orbit_reconst_20241009090000_20241009095959_20241010000530.nc||${L0_dir}/prefire_02_orbit_reconst_20241009100000_20241009105959_20241010000530.nc||${L0_dir}/prefire_02_orbit_reconst_20241009110000_20241009115959_20241010000530.nc||${L0_dir}/prefire_02_orbit_reconst_20241009120000_20241009125959_20241010000530.nc||${L0_dir}/prefire_02_orbit_reconst_20241009130000_20241009135959_20241010000530.nc||${L0_dir}/prefire_02_orbit_reconst_20241009140000_20241009145959_20241010000530.nc||${L0_dir}/prefire_02_orbit_reconst_20241009150000_20241009155959_20241010000544.nc||${L0_dir}/prefire_02_orbit_reconst_20241009160000_20241009165959_20241010000544.nc||${L0_dir}/prefire_02_orbit_reconst_20241009170000_20241009175959_20241010000544.nc||${L0_dir}/prefire_02_orbit_reconst_20241009180000_20241009185959_20241010120534.nc||${L0_dir}/prefire_02_orbit_reconst_20241009190000_20241009195959_20241010120534.nc||${L0_dir}/prefire_02_orbit_reconst_20241009200000_20241009205959_20241010120534.nc||${L0_dir}/prefire_02_orbit_reconst_20241009210000_20241009215959_20241010120534.nc";


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

#====== For certain types of development and testing ======
#SRF_DISAMBIG_STR="SRF_v12_2024-04-18";

#export SRF_DISAMBIG_STR;
#==========================================================

OUTPUT_DIR="${base_dir}/outputs/m3";

PROC_MODE=3;

  # * Only increment 'Rxx' when the resulting products will be DAAC-ingested
PRODUCT_FULLVER="R01_P00";
  # Special form ('R00_Syy') when processing simulated observations:
#PRODUCT_FULLVER="R00_S01";

export PACKAGE_TOP_DIR ANCILLARY_DATA_DIR INSTRUMENT_MODEL_DIR;
export OUTPUT_DIR PROC_MODE PRODUCT_FULLVER;

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

for cfg_str in ${L0_granule_cfg_str1} ${L0_granule_cfg_str2}
do
   GRANULE_START_ID_END=${cfg_str%:*};  # [0]

   tmp_str=${cfg_str##*:};  # [1]
   a=($(echo "$tmp_str" | tr ',' ' '))
   L0_PAYLOAD_FPATHS="${a[0]}";
   L0_BUS_FPATHS="${a[1]}";
   L0_ORBIT_FPATHS="${a[2]}";

   export GRANULE_START_ID_END L0_PAYLOAD_FPATHS L0_BUS_FPATHS;
   export L0_ORBIT_FPATHS;

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
