# ADIAB_WAVE demonstration case

We here provide a physical description of this experiment and additional details as to how to run this experiment within NEMO. This experiment is **created and tested** for NEMO **code at revision 15197**. 
A ipython notebook is also provided as a demonstration of possible analysis.

## Objectives

The purpose of this test case is to validate the implementation of the Generalized Lagrangian Mean equations for the coupling/forcing of NEMO with waves. This test case was first proposed by Ardhuin et al. (2008) and was succesively detailed by Bennis et al (2011).

## Physical description

The ADIAB_WAVE test case consists on a steady monochromatic wave shoaling on a slope from 4 to 6 meter depth. This test case does not take into account the wave breaking nor bottom friction nor wave-induced mixing and neglects the fluid viscosity. The domain is rectangular with a size of 780m x 50m and with an uniform horizontal grid spacing rn_dx=rn_dy=10m. The flow is confined to a channel with free-slip boundary condition at North and South and the bottom topography is uniform in the y-direction. The vertical grid is a terrain following sigma-coordinates with 120 layers. The domain is initialized at rest with a constant temperature and salinity.

In this testcase NEMO is forced offline by the wave model WaveWatch3 (WW3).
The initial characteristics of the wave field are a significant wave height Hs=1.02m, a wave period T=5.26s, and a wave direction theta=90° (propagation in the x direction).

The Breivik implementation used for the calculation of the Stokes drift vertical decay is based on the deep water approximation which is not valid in this test case. We added the Stokes drift vertical decay parameterization as in Michaud et al (2012), which is valid for intermediate/shalow-water, in the MY_SRC/sbcwave.F90. 

## Exemple of run

The Stokes drift vertical profil parameterization for intermediate water is activated by the option ln_STOKES_ADIAB= .TRUE. in the namusr_def section of the namelist.

~~~fortran
!-----------------------------------------------------------------------
&namusr_def    !   User defined :   ADIAB_WAVE configuration
!-----------------------------------------------------------------------
   !                        ! type of vertical coordinate
   ln_zco      = .false.    ! z-coordinate
   ln_zps      = .false.    ! z-partial-step coordinate
   ln_sco      = .true.     ! s-coordinate
   rn_dx       =   10.      ! horizontal resolution   [meters]
   rn_dy       =   10.
   rn_dz       =   0.05     ! vertical   resolution   [meters]
   ln_STOKES_ADIAB = .true.    ! Stokes Drift (Shallow/Intermediate water)
 /
~~~

The input_data folder contains the wave forcing as well as the East and West open boundary forcing.
The wave forcing file need to be uncompressed before running.

Run the executable :  
```
mpirun -np 1 ./nemo 
```

The ADIAB_GLM_grid_U.nc output file should be generated.

The jupyter-notebook adiabwave_notebook.ipynb is set up to generate the cross-section plots of the Stokes drift, the quasi-Eulerian velocity as well as the Lagrangian velocity after 900s of simulation. Results shown an irrotational flow. The Lagrangian velocity has a vertical shear that is entirely due to the Stokes drift and the quasi-Eulerian velocity is homogeneous over the water column.

## References
Ardhuin, F., Rascle, N., Belibassakis, K.A., 2008. Explicit wave-averaged primitive equations using a generalized Lagrangian mean. Ocean Model. 20, 35–60.

Bennis, A.C., Ardhuin, F., Dumas, F., 2011. On the coupling of wave and three-dimensional circulation models: Choice of theoretical framework, practical implementation and adiabatic tests. Ocean Model. 40, 260–272.

Michaud, H., Marsaleix, Leredde, P., Estournel, C., Bourrin, F., Lyard, F., Mayet, C., Ardhuin, F., 2012. Three-dimensional modelling of wave-induced current from the surfzone to the inner shelf. Ocean Sci. 8, 657–681.
