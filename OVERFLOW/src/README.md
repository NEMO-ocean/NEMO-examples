# this is where we store the actual NEMO code for the test case 
 - EXP00
 - MY_SRC
 
 **NOTE: for this test case MY_SRC is empty because all routine sare already in NEMO svn official repository**
 
 
 In EXP00 directory there is available some namelists :
These namelists have same blocks of namelist_ref with choice of:
- tracer advection scheme = **FCT2** or **FCT4**
 - FCT2 = COMPACT 2nd order on horizontal and vertical
 - FCT4 = COMPACT 4th order on horizontal and vertical
- form of the momentum advectionvector = **flux** or **form**
- momentum advection scheme = **cen2** or **ubs**
 - cen2 = 2nd order centered scheme
 - ubs = 3rd order UBS scheme
- ln_dynvor_ene = .false. !  enstrophy conserving scheme
- ln_dynvor_ens = .true.  !  energy conserving scheme
- ln_dynvor_mix = .false. !  mixed scheme
- ln_dynvor_een = .false. !  energy & enstrophy scheme