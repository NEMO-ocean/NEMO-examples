# SEAMOUNT demonstration case

## Objectives
The SEAMOUNT test case illustrates the impact of different horizontal pressure gradient (HPG) schemes on spurious velocities in the case of terrain following coordinates (or any non-geopotential coordinate).  This test case was first proposed by Beckmann and Haidvogel (1993) and has become the standard test case in the development of HPG schemes.

## Physical description
The SEAMOUNT test case comprises of an isolated SEAMOUNT defined as a radially symmetric Gaussian bump in the centre of an East-West periodic channel.  The domain is initialised at rest with a horizontally uniform vertical density profile, defined via the temperature field.  The default settings are for a 4050m SEAMOUNT in a 380km x 288km channel of depth 4500m, with uniform horizontal grid spacing ∆x=∆y=4000m.  The vertical grid is a terrain following s-coordinates with 10 layers equally spaced in the depth-normalised sigma dimension.  The analytical solution is for the ocean to remain at rest, any velocities are spurious.  The dimensions of the domain and the strength of the vertical density gradient are parameters which may be set in the namusr_def section of the namelist:
!-----------------------------------------------------------------------
&namusr_def    !    SEAMOUNT TEST CASE
!-----------------------------------------------------------------------
   rn_dx             = 4000.0       ! horizontal resolution  [m]
   rn_length         = 380.0        ! length of domain      [km]
   rn_width          = 288.0        ! width of domain       [km]
   rn_dz             = 450.0        ! vertical   resolution  [m]
   rn_initrho        = 0.1          ! Initial density perturbation magnitude
   rn_s              = 2.0          ! Burger number (to control background density profile)
   rn_bathy          = 4500.0       ! Max Depth
   rn_seamountheight = 4050.0       ! Seamount height
   rn_l              = 25000.0      ! Gaussian scale factor
   rn_f              = 0.0001       ! Coriolis parameter
   ln_exp_init       = .true.       ! Exponential decay density perturbation profile
   ln_linear_init    = .false.      ! Linearly decreasing density perturbation


## Experimental set-up
The intended experimental set up is for users to change the hpg scheme and compare the resulting hpg trend, and the evolution of the maximum spurious velocities over 10 days of spin-up.

!-----------------------------------------------------------------------
&namdyn_hpg    !   Hydrostatic pressure gradient option                 (default: NO selection)
!-----------------------------------------------------------------------
   ln_hpg_zco  = .false.   !  z-coordinate - full steps
   ln_hpg_zps  = .false.   !  z-coordinate - partial steps (interpolation)
   ln_hpg_sco  = .true.   !  s-coordinate (standard jacobian formulation)
   ln_hpg_isf  = .false.   !  s-coordinate (sco ) adapted to isf
   ln_hpg_djc  = .false.   !  s-coordinate (Density Jacobian with Cubic polynomial)
      ln_hpg_djc_vnh = .true.  !  hor.  bc type for djc scheme (T=von Neumann, F=linear extrapolation)
      ln_hpg_djc_vnv = .true.  !  vert. bc type for djc scheme (T=von Neumann, F=linear extrapolation)
   ln_hpg_prj  = .false.   !  s-coordinate (Pressure Jacobian scheme)

Users should also change the name of the run to match the hpg scheme selected:
!-----------------------------------------------------------------------
&namrun        !   parameters of the run
!-----------------------------------------------------------------------
   cn_exp      =  "SEAMOUNT_xxx"!  experience name

This should generate 4 .nc files called:

SEAMOUNT_xxx_1d_grid_[TUVW].nc

and a run.stat.nc file which contains the maximum spurious velocity per time-step.  The jupyter-notebook seamount_notebook.ipynb is set up to generate x-z contour plots of the absolute hpg trend over the SEAMOUNT for each scheme tested, as well as a difference to a reference.  It also generates a time series of the evolution of the spurious velocities.

## References
Beckmann, A., Haidvogel, D.B., 1993. Numerical simulation of flow around a tall isolated seamount. Part I: Problem formulation and accuracy. J. Phys. Oceanogr. 23, 1736–1753.



