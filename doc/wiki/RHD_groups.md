

This set of parameters, contained in the namelist block `&RT_GROUPS`, sets radiation group properties for RAMSES RHD runs. Note that the number of photon groups (`NGROUPS`) is a compilation parameter, to be set in the Makefile. The default settings are HI, HeI and HeII ionizing photon groups with ionisation cross sections and energies derived from a blackbody spectrum with an effective temperature of 10^5 Kelvin.

For detailed descriptions of the concepts described here, see the following papers:

* [[1] RAMSES-RT: radiation hydrodynamics in the cosmological context](http://arxiv.org/abs/1304.7126)
* [[2] A scheme for radiation pressure and photon diffusion with the M1 closure in RAMSES-RT](http://arxiv.org/abs/1411.6440)

 
| Variable name, syntax, default value | Fortran type  | Description       |
|:---------------------------- |:------------- |:------------------------- |
| `groupL0 = 13.60,24.59,54.42` |  `real array`       | Lower energy boundaries, in eV, of each photon group (see §2 in [1]). Used for calculating SED model emission from stellar particles.|
| `groupL1 = 24.59,54.42,0.0`   |  `real array`       | Upper energy boundaries, in eV, of each photon group. A value of 0.0 is used to represent infinity.|
| `group_egy = 18.85,35.079,65.666` |  `real array`       | Average photon energies (eV) for each group. These can either be set manually or left to RAMSES to derive from SED models  (with `SEDprops_update>0`).|
| `group_csn(1,:) = 3.0d-18,0.0,0.0`          `group_csn(2,:) = 5.7d-19,4.5d-18,0.0` `group_csn(3,:) = 7.9d-20,1.2d-18,1.1d-18` |  `real matrix`    | 2d matrix representing average ionisation cross sections (cm^2) between each group (first index) and species (second index). These can either be set manually or left to RAMSES to derive from SED models (with `SEDprops_update>0`). |
| `group_cse(1,:) = 2.8d-18,0.0,0.0`  `group_cse(2,:) = 5.0d-19,4.1d-18,0.0` `group_cse(3,:) = 7.4d-20,1.1d-18,1.0d-18` |  `real matrix`    | 2d matrix representing energy weighted average ionisation cross sections (cm^2) between each group (first index) and species (second index). These can either be set manually or left to RAMSES to derive from SED models (with `SEDprops_update>0`). |
| `spec2group = 1,2,3` |  `integer array`       | Determines, for each recombining species (HII, HeII, HeIII) which photon group the recombination photons are injected into. Note that recombination emission must be activated with `rt_otsa=.false.` in `$RT_PARAMS`.|
| **====================** | **====================** | **Radiation pressure and IR radiation parameters** | 
| `kappaAbs = 0.0` |  `real array`       | Dust absorption (Planck) opacity factor. The real opacity scales with the local metallicity and if `is_kIR_T=.true.`, the IR opacity also varies with the local gas temperature.|
| `kappaSc = 0.0` |  `real array`       | Scattering (Rosseland) opacity factor, only used for the IR photon group (which is the first group). The real opacity scales with the local metallicity and if `is_kIR_T=.true.`, it also varies with the local gas temperature.|