# PREFIRE_L1

Package (written primarily in MatLab and Python) to produce the PREFIRE 1A-RAD and 1B-RAD science data products. These products provide spacecraft observation geometry and geolocated, calibrated TIRS radiances and brightness temperatures.

This code is released under the terms of this [LICENSE](LICENSE).  The version of this package can be found in [VERSION.txt](VERSION.txt).

# Installation

## Requirements

Python version 3.8+ is required, along with the following third-party Python packages: numpy, netcdf4, scipy

The associated (Python-based) git repositories ['PREFIRE_tools'](https://github.com/UW-PREFIRE/PREFIRE_tools) and ['PREFIRE_PRD_GEN'](https://github.com/UW-PREFIRE/PREFIRE_PRD_GEN) are also required for the proper operation of this package.

## Python Environment Setup

It is recommended to install the above Python packages in a dedicated conda environment (or something similar).  An example of the packages (and their versions) used successfully by us can be found in [conda_env.list](conda_env.list).

For example, using conda (and specifying Python 3.10.x from the conda-forge channel):

```
conda create --name for_PREFIRE_L1 -c conda-forge python=3.10;
conda activate for_PREFIRE_L1;
conda install -c conda-forge numpy netcdf4 scipy;
```

The location of 'PREFIRE_PRD_GEN' and 'PREFIRE_tools' depends on the value of the user's PYTHONPATH and/or sys.path -- for example, one could simply add each of those git repositories' local root Python source code directory to PYTHONPATH. Operationally, however, this package uses symbolic links to those git repositories' local root Python source code directories (or full copies of the same) in the source/ directory.

## Environment Variables

### Each job (executing this science algorithm package) is configured via information contained within environment variables.

### To specify that numpy, scipy, et cetera used by this algorithm should not use more than one thread or process, the below environment variables are expected to be set:

```
MKL_NUM_THREADS=1
NUMEXPR_NUM_THREADS=1
OMP_NUM_THREADS=1
VECLIB_MAXIMUM_THREADS=1
OPENBLAS_NUM_THREADS=1
```

### Some environment variables are always required to be set (also see test/run*.sh or test/run*.ps1):

PACKAGE_TOP_DIR  :  the top-level directory (i.e., the one that contains dist/, test/, etc.) of this package

ANCILLARY_DATA_DIR  :  the package's ancillary data directory (should be an absolute path)

INSTRUMENT_MODEL_DIR  :  the instrument model directory, typically a subdirectory of the package's ancillary data directory (should be an absolute path)

OUTPUT_DIR  :  the directory in which all meaningful output will be written (should be an absolute path)

PRODUCT_FULLVER  :  the full product processing/revision version string (e.g., "R01_P02").  Operationally, we only increment 'Rxx' when the resulting products will be DAAC-ingested.

PROC_MODE  :  the processing mode (operationally-valid values: 3, 4; see below for more details)

### This science algorithm package has two different operational modes, specified by the current value of the PROC_MODE environment variable (see above).  Each mode has a set of mode-specific environment variables.  More information about the values of each of these can be found in the `test/run*.sh` scripts.

== _PROC_MODE = 3  (produce 1A-RAD product granule from prepped L0-payload, L0-bus, and L0-orbit files)_

GRANULE_START_ID_END  :  granule ID string and the start/end timestamps of the granule

L0_PAYLOAD_FPATHS  :  filepath(s) of the prepped L0-payload telemetry file(s) that contain this granule, 17 hours before it, and 1.6 hours (~1 orbit) after it -- for calibration purposes

L0_BUS_FPATHS  :  filepath(s) of the prepped L0-bus telemetry file(s) that contain this granule, 17 hours before it, and 1.6 hours (~1 orbit) after it -- for calibration purposes

L0_ORBIT_FPATHS  :  filepath(s) of the prepped L0-orbit reconstruction file(s) that contain this granule, 17 hours before it, and 1.6 hours (~1 orbit) after it -- for calibration purposes

SRF_DISAMBIG_STR  :  [optional, for certain types of development and testing] a disambiguation string specifying which SRF/NEdR file to use (out of a selection of files in dist/ancillary/SRF/)

== _PROC_MODE = 4  (produce all or a portion of a 1B-RAD product granule, given a 1A-RAD product granule file)_

L1A_RAD_FILE  :  filepath of the "source" 1A-RAD product granule

ATRACK_IDX_RANGE_0BI  :  coded frame (i.e., along-track segment) subset to process and output (for example, "ATRACK_IDXRANGE_0BASED_INCLUSIVE:2001:3100" => atrack dim indices, from 2001 through 3100; "ATRACK_IDXRANGE_0BASED_INCLUSIVE:0:END" => atrack dim indices, 0 through the last frame)

DEM_ROOT_DIR  :  path (ending with `/tiles`) where the files of the tiled Copernicus GLO-90 (global, 90 m) digital elevation map (DEM) are stored

# Running the test script(s)

## Obtain and unpack ancillary and test data

## Spectral Response Functions (SRFs)

Spectral Response Function data (v13_2024-09-15) for the two Thermal InfraRed Spectrometer (TIRS) instruments aboard the dual NASA PREFIRE mission CubeSats are needed for this software package.  A zip-archive file containing that data can be downloaded from [here](https://zenodo.org/records/16638853).

To install the downloaded SRF data:

(1) If needed, create a subdirectory within this software package (`dist/ancillary/SRF/`).

(2) Then extract (unzip) the two SRF files from the downloaded zip-archive file, and put them both in that `dist/ancillary/SRF/` directory.

### Copernicus GLO-90 DEM

To create full 1B-RAD data products, the full Copernicus GLO-90 (global, 90 m) tiled digital elevation map dataset (nearly 200 GB in size) is required to be downloaded and present in a directory on the local system (see the DEM_ROOT_DIR environment variable description above).

Information about obtaining this dataset can be found [here](https://dataspace.copernicus.eu/explore-data/data-collections/copernicus-contributing-missions/collections-description/COP-DEM).

### Prepare the test input and output directories:

`cd test;`

On Linux/UNIX systems, possibly create a useful symbolic link to the test input data (if needed):

`ln -s WHEREEVER_THE_DATA_IS/inputs inputs;`

Prepare the output directory (Linux/UNIX example):

`mkdir -p outputs;`

_OR_ perhaps something like

`ln -s /data/users/myuser/data-PREFIRE_L1/outputs outputs;`

## Run the L1 package

### A Linux/UNIX example

`cp run.sh my-run.sh;`

Edit `my-run.sh` as needed (e.g., change input file names)

`./my-run.sh`

The output file(s) will be in subdirectories of `test/outputs/` (e.g., `m4/`)

## _The creation of this code was supported by NASA, as part of the PREFIRE (Polar Radiant Energy in the Far-InfraRed Experiment) CubeSat mission._
