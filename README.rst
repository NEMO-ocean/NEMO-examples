**********************
Explore the test cases
**********************

.. todo::

   CANAL animated gif is missing

.. contents::
   :local:
   :depth: 1

Installation
============

Download
--------

| The complete and up-to-date set of test cases is available on
  :github:`NEMO test cases repository <NEMO-examples>`.
| Download it directly into the :file:`./tests` root directory with

.. code-block:: console

   $ git clone http://github.com/NEMO-ocean/NEMO-examples

Compilation
-----------

The compilation of the test cases is very similar to
the manner the reference configurations are compiled.
If you are not familiar on how to compile NEMO,
it is first recomended to read :doc:`the instructions <install>`.

| As the reference configurations are compiled with ``-r`` option,
  test cases can be compiled by the use of :file:`makenemo` with ``-a`` option.
| Here an example to compile a copy named WAD2 of the wetting and drying test case (WAD):

.. code-block:: console

   $ ./makenemo -n 'WAD2' -a 'WAD' -m 'my_arch' -j '4'

Run and analysis
----------------

There no requirement of specific input file for the test_cases presented here.
The XIOS xml input files and namelist are already setup correctly.
For detailed description and Jupyter notebook, the reader is directed on
the :github:`NEMO test cases repository <NEMO-examples>`

The description below is a brief advertisement of some test cases.

List of test cases
==================

ICE_AGRIF
---------

.. figure:: _static/ICE_AGRIF_UDIAG_43days_UM5.gif
   :width: 200px
   :align: left

   ..

| This test case illustrates the advection of an ice patch across
  an East/West and North/South periodic channel over a slab ocean (i.e. one ocean layer),
  and with an AGRIF zoom (1:3) in the center.
| The purpose of this configuration is to
  test the advection of the ice patch in and across the AGRIF boundary.
  One can either impose ice velocities or ice-atm.
  Stresses and let rheology define velocities (see :file:`README` for details)

VORTEX
------

.. figure:: _static/VORTEX_anim.gif
   :width: 200px
   :align: right

   ..

This test case illustrates the propagation of an anticyclonic eddy over a Beta plan and a flat bottom.
It is implemented here with an online refined subdomain (1:3) out of which the vortex propagates.
It serves as a benchmark for quantitative estimates of nesting errors as in :cite:`DEBREU2012`,
:cite:`PENVEN2006` or :cite:`SPALL1991`.

The animation (sea level anomaly in meters) illustrates with
two 1:2 successively nested grids how the vortex smoothly propagates out of the refined grids.

ISOMIP
------

.. figure:: _static/ISOMIP_moc.png
   :width: 200px
   :align: left

   ..

| The purpose of this test case is to evaluate the impact of various schemes and new development with
  the iceshelf cavities circulation and melt.
  This configuration served as initial assesment of the ice shelf module in :cite:`LOSCH2008` and
  :cite:`MATHIOT2017`.
  The default setup is the one described |ISOMIP|_.
| The figure (meridional overturning circulation) illustrates
  the circulation generated after 10000 days by the ice shelf melting (ice pump).

.. |ISOMIP| replace:: here

LOCK_EXCHANGE
-------------

.. figure:: _static/LOCK-FCT4_flux_ubs.gif
   :width: 200px
   :align: right

   ..

| The LOCK EXCHANGE experiment is a classical fluid dynamics experiment that has been adapted
  by :cite:`HAIDVOGEL1999` for testing advection schemes in ocean circulation models.
  It has been used by several authors including :cite:`BURCHARD2002` and :cite:`ILICAK2012`.
  The LOCK EXCHANGE experiment can in particular illustrate
  the impact of different choices of numerical schemes and/or subgrid closures on
  spurious interior mixing.
| Here the animation of the LOCK_EXCHANGE test case using
  the advection scheme FCT4 (forth order) for tracer and ubs for dynamics.

OVERFLOW
--------

.. figure:: _static/OVF-sco_FCT4_flux_cen-ahm1000.gif
   :width: 200px
   :align: left

   ..

| The OVERFLOW experiment illustrates the impact of different choices of numerical schemes and/or
  subgrid closures on spurious interior mixing close to bottom topography.
  The OVERFLOW experiment is adapted from the non-rotating overflow configuration described in
  :cite:`HAIDVOGEL1999` and further used by :cite:`ILICAK2012`.
  Here we can assess the behaviour of the second-order tracer advection scheme FCT2 and
  forth-order FCT4, z-coordinate and sigma coordinate (...).
| Here the animation of the OVERFLOW test case in sigma coordinate with
  the forth-order advection scheme FCT4.

WAD
---

.. figure:: _static/wad_testcase_7.gif
   :width: 200px
   :align: right

   ..

| A set of simple closed basin geometries for testing the Wetting and drying capabilities.
  Examples range from a closed channel with EW linear bottom slope to
  a parabolic EW channel with a Gaussian ridge.
| Here the animation of the test case 7.
  This test case is a simple linear slope with a mid-depth shelf with
  an open boundary forced with a sinusoidally varying ssh.
  This test case has been introduced to emulate a typical coastal application with
  a tidally forced open boundary with an adverse SSH gradient that,
  when released, creates a surge up the slope.
  The parameters are chosen such that
  the surge rises above sea-level before falling back and oscillating towards an equilibrium position.

CANAL
-----

.. figure:: _static/CANAL_image.gif
   :width: 200px
   :align: left

   ..

East-west periodic canal of variable size with several initial states and
associated geostrophic currents (zonal jets or vortex).

ICE_ADV2D
---------

| This test case illustrates the advection of an ice patch across
  an East/West and North/South periodic channel over a slab ocean (i.e. one ocean layer).
  The configuration is similar to ICE_AGRIF, except for the AGRIF zoom.
| The purpose of this configuration is to test the advection schemes available in the sea-ice code
  (for now, Prather and Ultimate-Macho from 1st to 5th order),
  especially the occurence of overshoots in ice thickness

ICE_ADV1D
---------

| This experiment is the classical :cite:`SCHAR1996` test case ,
  which has been used in :cite:`LIPSCOMB2004`, and in which very specific shapes of ice concentration,
  thickness and volume converge toward the center of a basin.
  Convergence is unidirectional (in x) while fields are homogeneous in y.
| The purpose of this configuration is to
  test the caracteristics of advection schemes available in the sea-ice code
  (for now, Prather and Ultimate-Macho from 1st to 5th order),
  especially the constitency between concentration, thickness and volume,
  and the preservation of initial shapes.

.. rubric:: References

.. bibliography:: tests.bib
   :all:
   :style: unsrt
   :labelprefix: T

ICE_RHEO
--------
| 

BENCH
-----
| Benchmark configuration. Allow to run any configuration (including ORCA type or BDY) with idealized grid
  and initial state so it does not need any input file other than the namelists.
  As usual, all configuration changes can be done through the namelist. 
  We provide 3 example of namelist_cfg to mimic ORCA1, OR025 or ORCA12 configurations.
  By default do not produce any output file. An extensive description of BENCH will be abailable in 
  Irrmann et al. 2021.

CPL_OASIS
---------
| This test case checks the OASIS interface in OCE/SBC, allowing to set up 
  a coupled configuration through OASIS. See CPL_OASIS/README.md for more information.

DIA_GPU
---------
| This is a demonstrator of diagnostic DIAHSB ported to GPU using CUDA Fortran. 
  Memory communications between host and device are asynchronous given the device has that capability. 
  This experiment is target for ORCA2_ICE_PISCES

TSUNAMI
---------
| just use dynspg_ts to simulate the propagation of an ssh anomaly (cosinus) in a box configuration
  with flat bottom and jpk=2.

DONUT
-----
| Donut shaped configuration to test MPI decomposition with bdy.

C1D_ASICS
---------
| 

DOME
----
| 

ICB
----
| ICB is a very idealized configuration used to test and debug the icb module.
  The configuration is box with a shallow shelf (40m) on the east and west part of the domain 
  with a deep central trough (> 100m).
  ICB are generating using the test capability of the icb model along a E-W line (this can easily be tuned).

STATION_ASF
-----------
| this demonstration test case can be used to perform a sanity test of the SBCBLK interface of
  NEMO.  It will test all the bulk-parameterization algorithms using an idealized
  forcing that includes a wide range of *SSX / surface atmospheric state*
  conditions to detect potential error / inconsistencies.  Both a short report and
  boolean output: *passed* or *failed* is provided as an output.

SWG
---
| Square bassin blown with an analytical wind. Vertical structure allows only one mode
  associated with reduced gravity to develop. This configuration is based on Adcroft & Marshall 1998.
  Also run with RK3 time stepping. 

ADIAB_WAVE
----------
| The purpose of this test case is to validate the implementation of the Generalized Lagrangian Mean equations for the coupling of NEMO with waves. This test case was first proposed by Ardhuin et al. (2008) and was successively detailed by Bennis et al (2011).
