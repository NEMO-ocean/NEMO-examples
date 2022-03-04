# Coupling with OASIS test case
The CPL_OASIS test case allows to set up and check a basic coupling of NEMOto a simple TOYATM fake tmposhere through the OASIS coupler. A very limited number of fields are exchanged between NEMO and the TOYATM. The tests checks that the fields are indeed exchanged through OASIS and that the ATSSTSST field of sea surface temparture received by the TOYATM makes sense. If the test is sucessful, it states that the set up of NEMO-OASIS interface in the NEMO SBC module is working fine. 
<br>
We here provide a description of details of this experiment so as as how to run it and check test is sucessful. This experiment is **created and tested** for NEMO ** revision 12573 (to be replace by NEMO release 4.2 by end 2020)**. 

## Objectives
This test case enables the OASIS interface in NEMO (in the OCE/SBC module). A few fields are sent and received by NEMO and by TOYATM (the simplified "atmosphere"). The success of this test (see below ** Verification**) indicates that 
* the OASIS interface in NEMO is functionnal (some fields are sent and received)
* The sea surface temperature received by the TOYATM makes sense

This test case can be seen as a template to set up a coupling between NEMO and an atmospheric model through OASIS.

## Detailed description

This test case is a set up of NEMO (dynamics, sea-ice and biogeochemistry) on a global 2° grid (as in ithe ORCA2_ICE_PISCES reference configuration), except NEMO is here coupled to a "toyatm" through OASIS.

NEMO is running 160 time-steps (10 days). The coupling is done at each timestep (nn_fsbc=1 in namelist_cfg). The fields exchanged with the toyatm are defined in the &namsbc_cpl namelist in namelist_cfg file.


This test case requires:
* NEMO (no mofication of source code, from rev 12573 or higher (e.g. compatible with NEMO reease 4.2
* OASIS (need to be downloaded and compiled, with correct paths set in your arch file for NEMO
* TOYATM (the simple toy in place of an atmospheric model) (located in NEMO/tools directory, need to created executable using tools/maketools command)

This tests/CPL_OASIS directory contains all the need files:
* cpp_CPL_OASIS.fcm defining the active cpp keys for NEMO
* EXPREF directory containing all the input files: namelists, xml files, a template job tu run the test case, and a script to check the results and produce the report
  	 * More specifically, the fields exchanged by NEMO through OASIS ae defined in the namelist_cfg input file (see &namsbc_cpl variables)


## Building the CPL_OASIS test case
* Download and compile OASIS
* Build the NEMO executable for this CPL_OASIS test case. First you need to add the correct OASIS library path in your arch file in the %OASIS_HOME variable. Then, in your local NEMO root directory:
``` 
./makenemo -a CPL_OASIS -n MYCPL_OASIS -m "your arch file"
```
This makenemo command will create the test case in cfgs/MYCPL_OASIS
* Build the TOYATM executable
```
cd tools 
./maketools -n TOYATM -m "your arch file"
```

## Running the test case 
```
cd tests/MYCPL_OASIS/EXP00
cp ../../CPL_OASIS/job_run_CPL_TESTCASE .
Adapt the job_run_CPL_TESTCASE to your target computer and run it
```
In this directory the job_run_cpl_testcase contains all the steps to run the testcase. These steps are commented in the file.
** After adapting the headers for your batch system and the paths for the files **, run this script through the batch system of your target computer.

## Verification and validation of the test case
The script gen_report.sh located in the CPL_OASIS directory allows to check if the run came to a sucessful ending:
```
cd MYCPL_OASIS/EXP00
cp ../../CPL_OASIS/gen_report.sh .
./gen_report.sh
```
If the report is successful, a final check sould be done by visualising the ATSSTSST_toyatm_01.nc (using ncview or any other visualiser for NETCDF files) and comparing it to the reference ref_ATSSTSST_last_time_step.jpg image in the CPL_OASIS directory: the two visualisations must look alike.

reference ref_ATSSTSST_last_time_step.jpg
.. image:: ref_ATSSTSST_last_time_step.jpg


