# in this directory we store the actual NEMO code for the test case 
 - EXP00
 - MY_SRC
 
 **NOTE:** for this test case, at **revision 8097** MY_SRC is empty because all routines are already in NEMO svn official repository.
 
 
 In EXP00 directory, of this test case, at the revision 8097 of NEMO some namelists are availables : <br>
 
<span style="color:red"> better : You can see NEMO Reference Manual : http:///.....  for explication of the choice of namelist done in this directory?????</style>

Here there is a legend of some possible choices : 

- tracer advection scheme = **FCT2** or **FCT4**
 - FCT2 = COMPACT 2nd order on horizontal and vertical
 - FCT4 = COMPACT 4th order on horizontal and vertical
- form of the momentum advectionvector = **flux** or **form**
- momentum advection scheme = **cen2** or **ubs**
 - cen2 = 2nd order centered scheme
 - ubs = 3rd order UBS scheme
- ln\_dynvor\_ene  : enstrophy conserving scheme
- ln\_dynvor\_ens  : energy conserving scheme
- ln\_dynvor\_mix  : mixed scheme
- ln\_dynvor\_een  : energy & enstrophy scheme