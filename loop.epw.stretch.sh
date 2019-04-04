#!/bin/bash

varname="wspd"
units="ratio"

qsub -N "acc.stretch" -v gcm="ACCESS1-0",varname=$varname,units=$units run.stretch.morph.fac.pbs
qsub -N "can.stretch" -v gcm="CanESM2",varname=$varname,units=$units run.stretch.morph.fac.pbs
qsub -N "cnr.stretch" -v gcm="CNRM-CM5",varname=$varname,units=$units run.stretch.morph.fac.pbs
qsub -N "csi.stretch" -v gcm="CSIRO-Mk3-6-0",varname=$varname,units=$units run.stretch.morph.fac.pbs
qsub -N "gfd.stretch" -v gcm="GFDL-ESM2G",varname=$varname,units=$units run.stretch.morph.fac.pbs
qsub -N "hcc.stretch" -v gcm="HadGEM2-CC",varname=$varname,units=$units run.stretch.morph.fac.pbs
qsub -N "hes.stretch" -v gcm="HadGEM2-ES",varname=$varname,units=$units run.stretch.morph.fac.pbs
qsub -N "inm.stretch" -v gcm="inmcm4",varname=$varname,units=$units run.stretch.morph.fac.pbs
qsub -N "mir.stretch" -v gcm="MIROC5",varname=$varname,units=$units run.stretch.morph.fac.pbs
qsub -N "mri.stretch" -v gcm="MRI-CGCM3",varname=$varname,units=$units run.stretch.morph.fac.pbs
