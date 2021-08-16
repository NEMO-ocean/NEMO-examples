# ADIAB_WAVE demonstration case

## Objectives

The ADIAB_WAVE test case is used to validate the implementation of the Generalized Lagrangian Mean equation in case of wave-current interaction. This test case was first proposed by Ardhuin et al. (2008) and is commonly used to test the implementation of the wave current interaction in open ocean.

## Physical description

The ADIAB_WAVE test case consists on a steady monochromatic wave shoaling on a slope from 4 to 6 meter depth. The analytical solution is calculated for an inviscid fluid. This test case does not take into account the wave breaking nor bottom friction nor wave-induced mixing. The domain is rectangular with a size of 780m x 50m and with an uniform horizontal grid spacing rn_dx=rn_dy=10m. The flow is confined to a channel with free-slip boundary condition at North and South and the bottom topography is uniform in the y-direction. The vertical grid is a terrain following sigma-coordinates with 120 layers. The domain is initialized at rest with a constant temperature and salinity.

In this testcase NEMO is forced offline by WaveWatch3 (WW3).
The initial characteristics of the wave field are a significant wave height Hs=1.02m, a wave period T=5.26s, and a wave direction theta=90°.


The Breivik implementation used for the calculation of the Stokes drift decay is based on the deep water approximation. For this test case an additional Stokes Drift profile parameterization is added in MY_SRC/sbcwave.F90. 
The formulation for the Stokes drift decay is the one from Michaud et al, 2012 for intermediate/shallow water and is activated by by ln_STOKES_ADIAB= .TRUE. in the namusr_def section of the namelist

Two .nc files called ADIAB_GLM_grid_[UV].nc should be generated.

The jupyter-notebook adiabwave_notebook.ipynb is set up to generate the colormap plots of the Stokes drift profile, the quasi-Eulrian velocity as well as the Lagrangian velocity after 900s of simulation. In absence of dissipation, and given proper lateral boundary conditions the flow in wave shoaling over a bottom slope has to be irrotational. The reference solution exhibit a vertical shear that is entirely due to the Stokes drift and the quasi-Eulerian velocity is homogeneous over the water column.

## References
Ardhuin, F., Rascle, N., Belibassakis, K.A., 2008. Explicit wave-averaged primitive equations using a generalized Lagrangian mean. Ocean Model. 20, 35–60.

Bennis, A.C. , Ardhuin, F., Dumas, F., 2011. On the coupling of wave and three-dimensional circulation models: Choice of theoretical framework, practical implementation and adiabatic tests. Ocean Model. 40, 260–272.
