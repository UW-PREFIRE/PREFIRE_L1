#!/usr/bin/env bash

## IMPORTANT: Only run this script from the directory it resides in, i.e. with
##             ./run.sh    OR    bash run.sh

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

#set -ve;  # Exit on the first error, and print out commands as we execute them
set -e;  # Exit on the first error

# Determine the absolute path of the current working directory:
#  (this is typically the package dist/ directory)
readonly base_dir="$(absfpath ".")";

this_top_dir="$(absfpath "${base_dir}/..")";

# Change dir to the 'test/' directory:
test_path="${this_top_dir}/test";
cd "$test_path";

# Directory (or a symlink of the same name) 'test/inputs/' must already exist
#  (and contain raw L0 telemetry files) for this package:
inputs_path="$test_path"/inputs;

# Remove old L1A output files, create new ones, and then copy them to inputs/:
mnk="m3";
outputs_path="$test_path"/outputs/$mnk;
rm -rf "$outputs_path";
mkdir -p "$outputs_path";
./run_$mnk.sh;
cp "$outputs_path"/raw-PREFIRE_SAT1_1A-RAD_R01_P00_20241017183356_02035.nc "$inputs_path";
cp "$outputs_path"/raw-PREFIRE_SAT2_1A-RAD_R01_P00_20241009175747_02077.nc "$inputs_path";

# Remove old L1B output files and create new ones:
mnk="m4";
outputs_path="$test_path"/outputs/$mnk;
rm -rf "$outputs_path";
mkdir -p "$outputs_path";
./run_$mnk.sh;
