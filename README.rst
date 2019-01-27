**********************
Explore the test cases
**********************

.. contents::
	:local:

List
====

  The description below is a brief description of the test cases available in NEMO. 

ICE_AGRIF
---------
  
  This test case illustrates the advection of an ice patch across an East/West and North/South periodic channel
  over a slab ocean (i.e. one ocean layer), and with an AGRIF zoom (1:3) in the center
  The purpose of this configuration is to test the advection of the ice patch in  
  and across the AGRIF boundary
  One can either impose ice velocities or ice-atm. stresses and let rheology define velocities
  (see README for details)

  .. image:: _static/ICE_AGRIF_UDIAG_43days_UM5.gif

VORTEX
------
  
  This test case illustrates the propagation of an anticyclonic eddy over a Beta plan and a flat bottom.
  It is implemented here with an online refined subdomain (1:3) out of which the vortex propagates.
  It serves as a benchmark for quantitative estimates of nesting errors as in Debreu et al. (2012) :cite:`DEBREU2012`,
  Penven et al. (2006) :cite:`PENVEN2006` or Spall and Holland (1991) :cite:`SPALL1991`.
  
  The animation below (sea level anomaly in meters) illustrates with two 1:2 successively nested grids how
  the vortex smoothly propagates out of the refined grids.
  
  .. image:: _static/VORTEX_anim.gif

ISOMIP
------

  The purpose of this test case is to evaluate the impact of various schemes and new development with the iceshelf cavities circulation and melt.
  This configuration served as initial assesment of the ice shelf module in Losh et al. (2008) :cite:`LOSCH2008` and Mathiot et al. (2017) :cite:`MATHIOT2017`. 
  The default setup is the one described `here <http://staff.acecrc.org.au/~bkgalton/ISOMIP/test_cavities.pdf>`_.
  
  The figure below (meridional overturning circulation) illustrates the circulation generated after 10000 days by the ice shelf melting (ice pump).

  .. image:: _static/ISOMIP_moc.png

LOCK_EXCHANGE
-------------

  The LOCK EXCHANGE experiment is a classical fluid dynamics experiment that has been adapted
  by Haidvogel and Beckmann (1999) :cite:`HAIDVOGEL1999` for testing advection schemes in ocean circulation models.
  It has been used by several authors including Burchard and Bolding (2002) :cite:`BURCHARD2002` and Ilicak et al. (2012) :cite:`ILICAK2012`.
  The LOCK EXCHANGE experiment can in particular illustrate the impact of different choices of numerical schemes 
  and/or subgrid closures on spurious interior mixing.

  Below the animation of the LOCK_EXCHANGE test case using the advection scheme FCT4 (forth order) for tracer and ubs for dynamics.

  .. image:: _static/LOCK-FCT4_flux_ubs.gif

OVERFLOW
--------

  The OVERFLOW experiment illustrates the impact of different choices of numerical schemes 
  and/or subgrid closures on spurious interior mixing close to bottom topography. 
  The OVERFLOW experiment is adapted from the non-rotating overflow configuration described 
  in Haidvogel and Beckmann (1999) :cite:`HAIDVOGEL1999` and further used by Ilicak et al. (2012) :cite:`ILICAK2012`.
  Here we can assess the behaviour of the second-order tracer advection scheme FCT2 and fortht-order FCT4, z-coordinate and sigma coordinate (...).

  Below the animation of the OVERFLOW test case in sigma coordinate with the forth-order advection scheme FCT4.

  .. image:: _static/OVF-sco_FCT4_flux_cen-ahm1000.gif

WAD
---

  A set of simple closed basin geometries for testing the Wetting and drying capabilities. 
  Examples range from a closed channel with EW linear bottom slope to a parabolic EW channel with a Gaussian ridge.
  
  Below the animation of the test case 7. This test case is a simple linear slope with a mid-depth shelf with an open boundary forced with a sinusoidally varying ssh.
  This test case has been introduced to emulate a typical coastal application with a tidally forced open boundary with an adverse SSH gradient that, when released, creates a surge up the slope.
  The parameters are chosen such that the surge rises above sea-level before falling back and oscillating towards an equilibrium position

  .. image:: _static/wad_testcase_7.gif

CANAL
-----

  East-west periodic canal of variable size with several initial states and associated geostrophic currents (zonal jets or vortex).

  .. image::_static/CANAL_image.gif

ICE_ADV2D
---------
  
  This test case illustrates the advection of an ice patch across an East/West and North/South periodic channel
  over a slab ocean (i.e. one ocean layer).
  The configuration is similar to ICE_AGRIF, except for the AGRIF zoom.
  The purpose of this configuration is to test the advection schemes available in the sea-ice code
  (for now, Prather and Ultimate-Macho from 1st to 5th order),
  especially the occurence of overshoots in ice thickness
  

ICE_ADV1D
---------
  
  This experiment is the classical Schar & Smolarkiewicz (1996) test case :cite:`SCHAR1996`,
  which has been used in :cite:`LIPSCOMB2004`,
  and in which very specific shapes of ice concentration, thickness and volume converge toward the center of a basin.
  Convergence is unidirectional (in x) while fields are homogeneous in y.
  The purpose of this configuration is to test the caracteristics of advection schemes available in the sea-ice code
  (for now, Prather and Ultimate-Macho from 1st to 5th order),
  especially the constitency between concentration, thickness and volume, and the preservation of initial shapes.
  
  

  
Compile test cases
==================

The compilation of the test cases is very similar to the manner the reference configurations are compiled.
If you are not familiar on how to compile NEMO, it is first recomended to read :doc:`the instructions <install>`

| In the same manner as the ref. cfg are compiled with '-r' option, test cases can be compile by the use of makenemo with '-a' option.

| Here an example to compile a copy named WAD2 of the wetting and drying test case (WAD) on the macport_osx architecture on 4 cores:

.. code-block:: console
 
	$ ./makenemo -n WAD2 -a WAD -m macport_osx -j 4

Run and analyse the test cases
==============================

There no requirement of specific input file for the test_cases presented here. The XIOS xml input files and namelist are already setup correctly. 
For detailed description and Jupyter notebook, the reader is directed on
the `NEMO test cases repository <http://github.com/NEMO-ocean/NEMO-examples>`_

References
==========

.. bibliography:: test_cases.bib
	:all:
	:style: unsrt
	:labelprefix: T
