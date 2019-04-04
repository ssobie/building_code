#!/bin/bash

qsub -N "acc.dwpt" -v gcm="ACCESS1-0" run.dwpt.morph.fac.pbs
qsub -N "can.dwpt" -v gcm="CanESM2" run.dwpt.morph.fac.pbs
qsub -N "cnr.dwpt" -v gcm="CNRM-CM5" run.dwpt.morph.fac.pbs
qsub -N "csi.dwpt" -v gcm="CSIRO-Mk3-6-0" run.dwpt.morph.fac.pbs
qsub -N "gfd.dwpt" -v gcm="GFDL-ESM2G" run.dwpt.morph.fac.pbs
qsub -N "hcc.dwpt" -v gcm="HadGEM2-CC" run.dwpt.morph.fac.pbs
qsub -N "hes.dwpt" -v gcm="HadGEM2-ES" run.dwpt.morph.fac.pbs
qsub -N "inm.dwpt" -v gcm="inmcm4" run.dwpt.morph.fac.pbs
qsub -N "mir.dwpt" -v gcm="MIROC5" run.dwpt.morph.fac.pbs
qsub -N "mri.dwpt" -v gcm="MRI-CGCM3" run.dwpt.morph.fac.pbs
