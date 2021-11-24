# ADIAB_WAVE demonstration case

We here provide a physical description of this experiment and additional details as to how to run this experiment within NEMO. This experiment is **created and tested** for NEMO **code at revision 15197**. 
An ipython notebook is also provided as a demonstration of possible analysis.

## Objectives

The purpose of this test case is to validate the implementation of the Generalized Lagrangian Mean equations for the coupling of NEMO with waves. This test case was first proposed by Ardhuin et al. (2008) and was successively detailed by Bennis et al (2011).

## Physical description

The ADIAB_WAVE test case consists on a steady monochromatic wave shoaling on a slope from 4 to 6-meter depth. It does not take into account the wave breaking nor bottom friction nor wave-induced mixing. The domain is rectangular with a size of 780m x 50m and with a uniform horizontal grid spacing rn_dx=rn_dy=10m. The flow is confined to a channel with a free-slip boundary condition at North and South and the bottom topography is uniform in the y-direction. The vertical grid is a terrain-following sigma-coordinates with 120 layers. The domain is initialized at rest with a constant temperature and salinity.

In a sake of simplicity, for this testcase, NEMO is forced offline by the wave model WaveWatch3 (WW3).
The initial characteristics of the wave field are a significant wave height Hs=1.02m, a wave period T=5.26s, and a wave direction theta=90° (propagation in the x-direction).

The calculation of the Stokes drift vertical decay in NEMO is currently based on the deep water approximation that is not verified in this test case. The MY_SRC/sbcwave.F90 includes the Stokes drift vertical decay parameterization for intermediate/shallow-water (Michaud et al, 2012).  

## Exemple of run

The calculation of the Stokes drift for intermediate water is activated by the option ln_STOKES_ADIAB= .TRUE. in the namusr_def section of the namelist.

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
The wave forcing file (sdw_WW3_ADIABWAVE_y0001m01d01.nc.tar.gz) need to be uncompressed before running.

Run the executable :  
```
mpirun -np 1 ./nemo 
```

The ADIAB_GLM_grid_U.nc output file should be generated.

The jupyter-notebook adiabwave_notebook.ipynb is set up to generate the vertical cross-section plots of the Stokes drift, the quasi-Eulerian velocity as well as the Lagrangian velocity after 900s of simulation. The quasi-Eulerian velocity is irrotational and homogeneous over the water column. The Lagrangian velocity vertical shear is entirely due to the Stokes drift.

## References
Ardhuin, F., Rascle, N., Belibassakis, K.A., 2008. Explicit wave-averaged primitive equations using a generalized Lagrangian mean. Ocean Model. 20, 35–60.

Bennis, A.C., Ardhuin, F., Dumas, F., 2011. On the coupling of wave and three-dimensional circulation models: Choice of theoretical framework, practical implementation and adiabatic tests. Ocean Model. 40, 260–272.

Michaud, H., Marsaleix, Leredde, P., Estournel, C., Bourrin, F., Lyard, F., Mayet, C., Ardhuin, F., 2012. Three-dimensional modelling of wave-induced current from the surfzone to the inner shelf. Ocean Sci. 8, 657–681.
