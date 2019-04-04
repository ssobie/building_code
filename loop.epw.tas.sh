#!/bin/bash

qsub -N "acc.tas" -v gcm="ACCESS1-0" run.tas.morph.fac.pbs
qsub -N "can.tas" -v gcm="CanESM2" run.tas.morph.fac.pbs
qsub -N "cnr.tas" -v gcm="CNRM-CM5" run.tas.morph.fac.pbs
qsub -N "csi.tas" -v gcm="CSIRO-Mk3-6-0" run.tas.morph.fac.pbs
qsub -N "gfd.tas" -v gcm="GFDL-ESM2G" run.tas.morph.fac.pbs
qsub -N "hcc.tas" -v gcm="HadGEM2-CC" run.tas.morph.fac.pbs
qsub -N "hes.tas" -v gcm="HadGEM2-ES" run.tas.morph.fac.pbs
qsub -N "inm.tas" -v gcm="inmcm4" run.tas.morph.fac.pbs
qsub -N "mir.tas" -v gcm="MIROC5" run.tas.morph.fac.pbs
qsub -N "mri.tas" -v gcm="MRI-CGCM3" run.tas.morph.fac.pbs
