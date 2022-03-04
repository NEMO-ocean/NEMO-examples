# *Station Air-Sea Fluxes* demonstration case

Last successful test done with NEMOGCM branch: `NEMO/branches/2020/dev_r13648_ASINTER-04_laurent_bulk_ice`, revision `13806`

Author: Laurent Brodeau, 2020

## Objectives

`STATION_ASF` is a demonstration test-case that mimics a (static) in-situ station (buoy, platform) dedicated to the estimation of surface air-sea fluxes by means of *widely-measured* (bulk) meteorological surface parameters.

`STATION_ASF` has been constructed by merging the *single column* and the *standalone surface module* configurations of NEMO. In short, it can be defined as "SAS meets C1D". As such, the spatial domain of `STATION_ASF` is punctual (1D, well actually 3 x 3 as in C1D).

`STATION_ASF` is therefore a versatile tool, and extremely lightweight in terms of computing requirements, to test the different bulk algorithms and cool-skin/warm-layer parameterization options included in NEMO.

As input `STATION_ASF` will require the traditional *bulk* sea surface parameters:
- Bulk sea surface temperature (SST) at _z<sub>SST</sub>_ meters below the surface
- Surface current vector
- Sea surface salinity

as well as the usual surface atmospheric state:
- air temperature at _z<sub>t</sub>_ meters above the surface
- air humidity  at _z<sub>t</sub>_ meters above the surface (specific humidity or relative humidity or dew-point temperature)
- wind speed vector at _z<sub>u</sub>_ meters above the surface
- Sea level atmospheric pressure (SLP)
- Downwelling solar radiation
- Downwelling longwave radiation


### Example of diagnostics with `STATION_ASF`

(Generated with script `./EXPREF/plot_station_asf_simple.py`)

![plot](https://github.com/NEMO-ocean/NEMO-examples/blob/master/STATION_ASF/figs/01_temperatures_ECMWF.svg)

![plot](https://github.com/NEMO-ocean/NEMO-examples/blob/master/STATION_ASF/figs/Cd.svg)

![plot](https://github.com/NEMO-ocean/NEMO-examples/blob/master/STATION_ASF/figs/dT_skin.svg)

![plot](https://github.com/NEMO-ocean/NEMO-examples/blob/master/STATION_ASF/figs/Qlat.svg)


## Physical description

### Important namelist parameters specific to `STATION_ASF`

* ```rn_dept1@namusr_def:``` depth (m) at which the prescribed SST is taken (*i.e.* depth of first T-point); important due to impact on warm-layer estimate, the deeper, the more pronounced!

* ```rn_lat1d,rn_lon1d@namc1d:``` fixed coordinates of the location of the station (buoy, platform, etc).

* ```namsbc_blk:``` to be filled carefully, just as for "C1D", the prescribed surface ATMOSPHERIC state (files) are time series of shape 3x3 in space

* ```namsbc_sas:``` to be filled carefully, just as for "C1D", the prescribed surface OCEAN state (files) are time series of shape 3x3 in space



## Testing the sanity of the SBCBLK (bulk atmospheric forcing) interface of NEMO

`STATION_ASF` can be used to perform a sanity test of the SBCBLK interface of
NEMO.  It will test all the bulk-parameterization algorithms using an idealized
forcing that includes a wide range of *SSX / surface atmospheric state*
conditions to detect potential error / inconsistencies.  Both a short report and
boolean output: *passed* or *failed* is provided as an output.


### How is the validation test performed ?

Turbulent fluxes and transfer coefficients computed by NEMO with all different bulk parameterizations (via `STATION_ASF`), in all these idealized possible air-sea meteorological states, must be checked against a *supposedly trustable* reference solution, plus/minus a tolerance range not to deviate from (to allow the deviation from the reference inherent to each tested bulk parameterization). This reference solution, and this tolerance, are constructed independently and outside of NEMO. They are based on the mean, and twice the standard deviation, respectively, across all the members (1 member = 1 bulk parameterizations, N=5). This reference solution and tolerance range, just as the idealized input atmospheric and sea-surface forcings, are provided in the same NetCDF file:
`STATION_ASF/input_data/input_output_VALIDATION_IDEALIZED.nc`


### Performing the sanity test

First compile the `STATION_ASF` test-case as follows (compile with xios-2.5 support → check your ARCH file):

```./makenemo -a STATION_ASF -m <your_arch> -n STATION_ASF2 -j 4```

(successful compilation generates ```tests/STATION_ASF2/BLD/bin/nemo.exe```

Move to: ```tests/STATION_ASF/SANITY_CHECK/```

There, in script ``sanity_check_SBCBLK.sh``, double check that the variables given a value in the first lines are consistent with your setup (it should).

Launch the script:

    $ ./sanity_check_SBCBLK.sh


The report is both interactively shown in the standard output and spawned in the report file:

- `SBCBLK.success` the test passed :)
- `SBCBLK.fail` the test failed :(

In case of failure, read the rapport to know which algorithm/computed variables did not pass the test... You should also re-run the script `sanity_check_SBCBLK.sh` with the "more" argument:

    $ ./sanity_check_SBCBLK.sh more

 This will generate more output such as figure of time series for the tested variables.

### Dependencies

In order for the python script `analyze_output.py` to work, you need `Python 3` with the following modules:
- NumPy
- netCDF4
- (Matplotlib)

### Example of the report for a successful test


    ############ FINAL REPORT ############

     ***** Algorithm "ECMWF" PASSED sanity check !!!

     ***** Algorithm "NCAR" PASSED sanity check !!!

     ***** Algorithm "COARE3p0" PASSED sanity check !!!

     ***** Algorithm "COARE3p6" PASSED sanity check !!!

     ***** Algorithm "ANDREAS" PASSED sanity check !!!

     Test performed on the following NEMO prognostic variables:
      ==> qsb_oce, qla_oce, qlw_oce, taum, Cd_oce, Ce_oce

       ####################################
       ###    TEST PASSED FOR SBCBLK !  ###
       ####################################



## Playing with `STATION_ASF`


### Input forcing files to test STATION ASF

By default, `STATION_ASF` comes with two forcing sets:

* The `PAPA` forcing is an "open-ocean-only" hourly forcing, that corresponds to a real off-shore platform, the PAPA station (North-East Pacific, 50°N, 145°W). Hence the data provided for the PAPA test is actual data measured by the PAPA station.

* The `ERA5_NorthGreenland` forcing simulates a *virtual station* as it is extracted from the ECMWF ERA5 reanalysis right north of Greenland (84°N,36°W). This forcing is also hourly and allows for `STATION_ASF` to be used with sea-ice support as this particular location is ice-covered year-round.


#### `PAPA` forcing

One full year (2018) of processed hourly data from the PAPA station (buoy) are found in the 3 following files (found in `tests/STATION_ASF/input_data`):
- ```Station_PAPA_50N-145W_atm_hourly_y2018.nc```  → contains hourly surface atmospheric state
- ```Station_PAPA_50N-145W_precip_daily_y2018.nc``` → contains daily precipitation
- ```Station_PAPA_50N-145W_oce_hourly_y2018.nc``` → contains hourly sea surface state

For station PAPA (50.1 N, 144.9 W), air temperature and humidity are measured at 2.5 m, the wind speed at 4 m, and the SST at 1 m below the surface, hence the following namelist parameters are given:

    ...
    &namusr_def
        rn_dept1 =    1. 
    ...
    &namc1d
        rn_lat1d =  50.1 
        rn_lon1d = 215.1
    ...
    &namsbc_blk
        ...
        rn_zqt = 2.5
        rn_zu  = 4.
        ...
        ln_humi_rlh = .true.  !  humidity in PAPA data is relative humidity  [%]
        ...

The namelists for the `PAPA` forcing are located into `EXPREF/PAPA/oce/`



#### `ERA5_NorthGreenland` forcing

The entire forcing is found in one single file: `ERA5_NorthGreenland_surface_84N_-36E_1h_y2018.nc`.
Since we are dealing with a reanalysis of the ECMWF:

    ...
    &namusr_def
        rn_dept1 = 1. ! we assume ECMWF uses a prescribed bulk SST product at the standard depth of 1m...
    ...
    &namc1d
        rn_lat1d =  84. 
        rn_lon1d = 324.
    ...
    &namsbc_blk
        ...
        rn_zqt =  2.
        rn_zu  = 10.
        ...
        ln_humi_dpt = .true.  !  humidity in ERA5 is dew-point temperature [K]
        ...

The namelists for the `ERA5_NorthGreenland` forcing are located into `ERA5/oce+ice/` or `ERA5/oce`, depending whether you want to use sea-ice support or not.



### Compiling and launching STATION_ASF simulations

First compile the test-case as follows (compile with xios-2.5 support → check your ARCH file):

    $ ./makenemo -a STATION_ASF -m <your_arch> -n STATION_ASF2 -j 4

Move to `tests/STATION_ASF/EXPREF/`, there, you can use the script ``launch_sasf.sh`` to launch `STATION_ASF` simulations. You need to adapt the following variable to your environment in the script:

- ```CONFIG_BLD``` : the name of the sub-directory of `tests/` in which `STATION_ASF` has been compiled (default: `CONFIG_BLD="STATION_ASF2"`)

- ```FORCING``` : this is where you pick the forcing you want to use, and wether to use sea-ice support or not (only available for `ERA5_NorthGreenland` forcing)

- ```PROD_DIR``` :  Directory in which the simulation will be installed and run

- ```MPI_LAUNCH```: this is where you can adapt how the `nemo.exe` executable is launched on your system (default: `MPI_LAUNCH="mpirun -n 1"`)


Now, it's time to launch simulations: 

    $ ./launch_sasf.sh

If everything goes according to plan, ``launch_sasf.sh`` should have generated the 4 following output files into `${PROD_DIR}/output` (example for PAPA forcing, 4 algos tested):

    STATION_ASF-ANDREAS_PAPA_1h_20180101_20181231_gridT.nc
    STATION_ASF-COARE3p6_PAPA_1h_20180101_20181231_gridT.nc
    STATION_ASF-ECMWF_PAPA_1h_20180101_20181231_gridT.nc
    STATION_ASF-NCAR_PAPA_1h_20180101_20181231_gridT.nc


#### Post-processing and figures

Use the `plot_station_asf_OCE.py` script (Python 3) to generate a multitude of figures that compares the *open ocean* transfer coefficients and fluxes between the 4 algorithms, in our case with the `PAPA` forcing:

    $ ./plot_station_asf_OCE.py -d <PROD_DIR> -f PAPA

Example of a comparison plot:
![plot](https://github.com/NEMO-ocean/NEMO-examples/blob/master/STATION_ASF/figs/Cd_PAPA_open-ocean.svg)

If you have used a forcing with sea-ice support, such as `ERA5_NorthGreenland`, you can do the same for *sea-ice* transfer coefficients and fluxes, which will compare *air-ice* bulk algorithms:

    $ ./plot_station_asf_ICE.py -d <PROD_DIR> -f ERA5_NorthGreenland

Example of a comparison plot:
![plot](https://github.com/NEMO-ocean/NEMO-examples/blob/master/STATION_ASF/figs/Cd_ERA5_NorthGreenland_sea-ice.svg)
Note: when sea-ice support is enabled, typically by using the following `FORCING` configuration in script `launch_sasf.sh`:

    FORCING="ERA5_NorthGreenland" ; i_sea_ice=1 ; SFORC="ERA5_NorthGreenland_surface_84N_-36E_1h"
    
the fact that `i_sea_ice=1` trigers the computation of *air-ice* fluxes with different bulk *air-ice* algorithms, over the remaining open-ocean fraction the same *air-sea* algorithm is used: `ECMWF`.


-----------------------

#### TO FINISH, BROKEN !!!

Then, you can fire the Jupyter notebook [`station_asf_notebook.ipynb`](https://github.com/NEMO-ocean/NEMO-examples/blob/master/STATION_ASF/notebook/station_asf_notebook.ipynb) found into the `notebook` directory! In which you should update `cprod_dir` to the same path as `PROD_DIR` of `launch_sasf.sh`...

---

*/Laurent, November 2020.*

