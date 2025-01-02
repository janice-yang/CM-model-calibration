#!/bin/sh

##########################################################################
#
#	Batch GA calibrations on data from one cell prep under one protocol sequence
#	Project CM-model-calibration (https://github.com/janice-yang/CM-model-calibration)
#
##########################################################################

## Set BSUB (LSF) template and template fill-in info
realdata=true
if $realdata; then
lsf_templ=batch_run_GA_single_exp.lsf
else
lsf_templ=batch_run_GA_single.lsf
fi
account=ACCOUNT_NAME
user_email=EMAIL
matlab_ver=MATLAB_VERSION
# alloc / premium (1.5x)/ expressalloc (2h max)/ low
queue=QUEUE_NAME
# Wall clock limit for job (hour) - 3h per protocol for 5min run
walltime=RUN_TIME

# Number of CPU cores required
ncores=NCORES_INT

## GA settings - shared
# Whether data should be normalized
normalization=true
# Whether to fit AP, CaT, or both ('APCaT')
datatype="'APCaT'"
# Whether CaT should be scaled when calculating fitness
scaleCaT=false
# Population initialization rng seeds
startseed=0
endseed=9

## GA settings - pseudodata only
# How much noise to add to pseudodata
sigmaAP=0
sigmaCaT=0
## GA settings - experimental data only
# File specifications - examples in comments
filename=FILENAME_STR
#filename="'DMG242PF_Fselected.xlsx'"
sheetnames=SHEETNAMES_CELLSTR
#sheetnames="{'100Na_1Hz','Dofetilide0nM_0.5Hz','Dofetilide1nM_0.5Hz'}"
# Number of beats in dataset
nbeats=NBEATS_LISTSTR
#nbeats="[11,6,6]"

## Path to .list files of cells and protocol sequences - change as needed
cell_list=cells.list
protocol_list=protocol_sequences.list

## Set up cell + protocol sequence
if $realdata; then
experiment=${filename//\'/}
experiment=${experiment//.xlsx/}
cell_ID=exp_${experiment}
else
cell_ID=$(cat $cell_list)
fi
protocol_sequence=$(cat $protocol_list)
# protocol_sequence is in format [#,#,#] - change to #-#-# for file name
sequence_name=${protocol_sequence//,/-}
sequence_name=${sequence_name//[/}
sequence_name=${sequence_name//]/}

##########################################################################

# Loop through rng seeds
for ((i=startseed; i<=endseed; i++))
do
# Random number generator seed
seed=$i

## Replace fillers in LSF template and create new job file
if [[ $normalization == "true" ]]; then
    jobname=$cell_ID.$sequence_name.norm.$seed.$datatype
else
    jobname=$cell_ID.$sequence_name.raw.$seed.$datatype
fi

if $realdata; then
sed "s|JBNAME|$jobname|g" $lsf_templ | \
sed "s|ACCOUNT|$account|g" | \
sed "s|USEREMAIL|$user_email|g" | \
sed "s|QUEUE|$queue|g" | \
sed "s|WALLTIME|$walltime|g" | \
sed "s|CORES|$ncores|g" | \
sed "s|VERSION|$matlab_ver|g" | \
sed "s|SEED|$seed|g" | \
sed "s|REALDATA|$realdata|g" | \
sed "s|FILENAME|$filename|g" | \
sed "s|PROTOCOLSEQ|$protocol_sequence|g" | \
sed "s|NORMALIZATION|$normalization|g" | \
sed "s|SCALECATBOOL|$scaleCaT|g" | \
sed "s|SHEETNAMES|$sheetnames|g" | \
sed "s|DATATYPE|$datatype|g" | \
sed "s|NBEATS|$nbeats|g" \
  > $jobname.lsf

else
sed "s|JBNAME|$jobname|g" $lsf_templ | \
sed "s|ACCOUNT|$account|g" | \
sed "s|USEREMAIL|$user_email|g" | \
sed "s|QUEUE|$queue|g" | \
sed "s|WALLTIME|$walltime|g" | \
sed "s|CORES|$ncores|g" | \
sed "s|VERSION|$matlab_ver|g" | \
sed "s|SEED|$seed|g" | \
sed "s|REALDATA|$realdata|g" | \
sed "s|CELLID|$cell_ID|g" | \
sed "s|PROTOCOLSEQ|$protocol_sequence|g" | \
sed "s|NORMALIZATION|$normalization|g" | \
sed "s|SCALECATBOOL|$scaleCaT|g" | \
sed "s|DATATYPE|$datatype|g" | \
sed "s|SIGMAAP|$sigmaAP|g" | \
sed "s|SIGMACAT|$sigmaCaT|g" \
  > $jobname.lsf
fi

#bsub < $jobname.lsf

done
