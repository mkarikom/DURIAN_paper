#!/bin/bash
# for baron

######################################################################################
# set project environment for R and julia
######################################################################################
nsleepsim=30 # amount of time to sleep after generating sc simulation (prevent pseudo from erroring upon creation)
nsleepfit=5 # amount of time to sleep in between steps that seem to miss key environment vars
nsleeploop=5 # how long to sleep between loop iterations
dsname=Baron
suffix=CPM
prefix=$dsname
SLURMACCT=qnie_lab

export EMDIAG=FALSE
# export SLURMPARTITION=highmem
export SLURMPARTITION=highmem
# export slurmtimelimit=3-00:00:00
export slurmtimelimit=1-00:00:00
export PROJECTDIR=/share/crsp/lab/cellfate/mkarikom/DURIAN_paper_clean
cd $PROJECTDIR
export BASEDIR=$PROJECTDIR/slurm
export OUTPUTMASTER=$BASEDIR/${dsname}/output.clusterMetrics.$SLURMPARTITION.$suffix
export NBULK=3
MEMP=16000M # memory in mb, try increasing if nodes are not avail
export NCPUS=$((NBULK+1))

export SOURCEPATH=$BASEDIR/${dsname}/durian_data

export JULIA_PROJECT=${PROJECTDIR} # make sure all workers can access the project enviroment
export JULIA_HOME=/opt/apps/julia/1.6.0/bin # make sure all workers can access the project enviroment
export R_HOME=/opt/apps/R/4.0.4/lib64/R # make sure JuliaCall/RCall can access R
export R_LIBS_USER=/data/homezvol2/mkarikom/R/x86_64-pc-linux-gnu-library/4.0 # norm() error
export JULIA_GR_PROVIDER=GR
export PYTHONPATH=/dfs6/pub/mkarikom/Python2.7_Pip_Packages

module purge
module load zlib
module load eigen
module load hdf5
module load gcc
module load lapack
module load OpenBLAS
module load mkl
module load julia/1.6.0
module load R/4.0.4
module load python/2.7.17 # needed for ursm, pypolyagamma
module load MATLAB/R2020b

export ETCLIB=$BASEDIR/scrabble_helper_functions/library_other_methods.R
export ALRALIB=$PROJECTDIR/ALRA/alra.R
export G2S3LIB=$PROJECTDIR/G2S3/run_G2S3.m
export CMFLIB=$PROJECTDIR/CMFImpute/analysis.m

# save all the enviroment stuff to the project dir
pip freeze > $PROJECTDIR/requirements.req
Rscript -e "sessionInfo()" >> $PROJECTDIR/session_info.req

######################################################################################
# durian/scrabble params
######################################################################################
export nEM=5 # the number of em iterations for durian and ursm
export ScrnIterOuter=10 # 10
export ScrnIterInner=10 # 10
export ScrnSDCIters=500000 # 500000
export DunIterOuter=10 # 10
export DunIterInner=10 # 10
export DunSDCIters=500000 # 500000
# export MCNITER=2500
export MCNITER=500

export ERR_IN_THRESH=1e-5 # 1e-5
export ERR_OUT_THRESH=1e-7 # 1e-7
export RUNOUTERSTATS=TRUE
export RUNSTABILITY=FALSE
export INITSCRABBLE=FALSE
export SUMMARIZEDECONV=TRUE
export USEIRLBA=FALSE # prevent instability for benchmarks, this will take longer

export DECONVGENETHRESH=-0.001
export SCRGENETHRESH=-0.001

SCPARAMS=( "1,1e-6,1e-4" "1e-2,1e-5,1e-5" )
DPARAMS=( "1,1e-6,1e-4" "1e-2,1e-5,1e-5" )

export durianEps=1e-3

export SUBSAMPLECPM=TRUE
export CELLSAMPLESIZE=1000
export SUBGENERATE=0.10

TypeList=( "" )
PREFIXTUPLES=( \
"BaronSC.DM.isletVST1K05;SegerstolpeBulk.DM.cpm" \
"BaronSC.H.isletVST1K05;SegerstolpeBulk.H.cpm" \
)

nsleepsim=60 # amount of time to sleep after generating sc simulation (prevent pseudo from erroring upon creation)
nsleepfit=2 # amount of time to sleep in between steps that seem to miss key environment vars
nsleeploop=0.05 # how long to sleep between loop iterations
export nsleepdatapath=1 # how long to sleep after creating pb data path
#################################################################


######################################################################################
# ursm params
######################################################################################
export number_of_cell_types=4 # for gupta is 11, for baron is 14
export burn_in_length=10
export gibbs_sample_number=10

for SUBSETCELLTYPES in "${TypeList[@]}"; do
        declare -a MULTISETDEPENDS=''
        for PREFIXTUPLE in "${PREFIXTUPLES[@]}"; do
                export SUBSETCELLTYPES
                export PREFIXTUPLE

                declare -a PSEUDODEPENDS=''
                declare -a NONURSMDEPENDS=''
                declare -a URSMDEPENDS2=''
                declare -a URSMDEPENDS3=''
                declare -a URSMDEPENDS4=''

                # export OUTBASEDIR=$OUTPUTMASTER/pref_$PREFIXTUPLE,dgene_$DECONVGENETHRESH,sgene_$SCRGENETHRESH,duEps_$durianEps,nEM_$nEM,sA_$ScrabbleAlpha,sB_$ScrabbleBeta,sG_$ScrabbleGamma,dA_$DurianAlpha,dB_$DurianBeta,dG_$DurianGamma,sub_$SUBSETCELLTYPES,suff_$suffix
                export OUTBASEDIR=$OUTPUTMASTER/pref_$PREFIXTUPLE,sub_$SUBSETCELLTYPES,suff_$suffix
                export DATAPATH=${OUTBASEDIR}/output_fit

                ######################################################################################
                # fit non-ursm imputation methods
                ######################################################################################

                SBATCHSUB=$BASEDIR/application_scripts/run_durian.sub
                IMPUTE_METHODS=( DrImpute dropout G2S3 CMFImpute )

                export JOBSCRIPT=$BASEDIR/application_scripts/run_imputation_methods.R

                for IMPUTE_METHOD in "${IMPUTE_METHODS[@]}"; do
                        export SIMREP=42
                        export IMPUTE_METHOD
                        echo running $IMPUTE_METHOD after $PSEUDODEPENDS
                        SBATCHJOBNAME=${PREFIXTUPLE}_${IMPUTE_METHOD}_$suffix
                        SBATCHOUTDIR=${OUTBASEDIR}/output_logs/pseudo_fit_$IMPUTE_METHOD/out
                        SBATCHERRDIR=${OUTBASEDIR}/output_logs/pseudo_fit_$IMPUTE_METHOD/err
                        mkdir -p $SBATCHERRDIR
                        mkdir -p $SBATCHOUTDIR
                        # collect the job ids `sbatchid` in an array
                        sbatchid=$(sbatch \
                        --account=$SLURMACCT \
                        --partition=$SLURMPARTITION \
                        --cpus-per-task=$NCPUS \
                        --time=$slurmtimelimit \
                        --mem-per-cpu=$MEMP \
                        --job-name=$SBATCHJOBNAME \
                        --error=$SBATCHERRDIR/err_%x_%A.log \
                        --out=$SBATCHOUTDIR/out_%x_%A.log \
                        $SBATCHSUB | cut -f 4 -d' ')
                        NONURSMDEPENDS+=":${sbatchid}"
                        MULTISETDEPENDS+=":${sbatchid}"
                        sleep $nsleeploop
                done

                ######################################################################################
                # fit scrabble and mtscrabble
                ######################################################################################
                IMPUTE_METHODS=( SCRABBLE )

                for IMPUTE_METHOD in "${IMPUTE_METHODS[@]}"; do
                        for SCPARAM in "${SCPARAMS[@]}"; do
                                IFS="," read -r ScrabbleAlpha ScrabbleBeta ScrabbleGamma <<< "${SCPARAM}"
                                export ScrabbleAlpha
                                export ScrabbleBeta
                                export ScrabbleGamma

                                export IMPUTE_METHOD
                                export SIMREP=42
                                echo running $IMPUTE_METHOD after $PSEUDODEPENDS
                                SBATCHJOBNAME=${PREFIXTUPLE}_${IMPUTE_METHOD}_$suffix
                                SBATCHOUTDIR=${OUTBASEDIR}/output_logs/pseudo_fit_$IMPUTE_METHOD/out
                                SBATCHERRDIR=${OUTBASEDIR}/output_logs/pseudo_fit_$IMPUTE_METHOD/err
                                mkdir -p $SBATCHERRDIR
                                mkdir -p $SBATCHOUTDIR
                                # collect the job ids `sbatchid` in an array
                                sbatchid=$(sbatch \
                                --account=$SLURMACCT \
                                --partition=$SLURMPARTITION \
                                --cpus-per-task=$NCPUS \
                                --time=$slurmtimelimit \
                                --mem-per-cpu=$MEMP \
                                --job-name=$SBATCHJOBNAME \
                                --error=$SBATCHERRDIR/err_%x_%A_%a.log \
                                --out=$SBATCHOUTDIR/out_%x_%A_%a.log \
                                $SBATCHSUB | cut -f 4 -d' ')
                                NONURSMDEPENDS+=":${sbatchid}"
                                MULTISETDEPENDS+=":${sbatchid}"
                                sleep $nsleeploop
                        done
                done

                ######################################################################################
                # fit durian music
                ######################################################################################

                for DPARAM in "${DPARAMS[@]}"; do
                        IFS="," read -r DurianAlpha DurianBeta DurianGamma <<< "${DPARAM}"

                        export DurianAlpha
                        export DurianBeta
                        export DurianGamma

                        export IMPUTE_METHOD=DURIAN
                        export DECONVMETHOD=MuSiC
                        export SIMREP=42

                        SBATCHJOBNAME=${PREFIXTUPLE}_${IMPUTE_METHOD}.${DECONVMETHOD}_$suffix
                        SBATCHOUTDIR=${OUTBASEDIR}/output_logs/pseudo_fit_$IMPUTE_METHOD.$DECONVMETHOD/out
                        SBATCHERRDIR=${OUTBASEDIR}/output_logs/pseudo_fit_$IMPUTE_METHOD.$DECONVMETHOD/err
                        mkdir -p $SBATCHERRDIR
                        mkdir -p $SBATCHOUTDIR

                        export NJULIACORES=$((NCPUS-1)) # this should be leq the number of PBULKS 
                        SBATCHSUB=$BASEDIR/application_scripts/pseudo_array_task.sub
                        export JOBSCRIPT=$BASEDIR/application_scripts/run_imputation_methods.R
                        export nCoresAvail=$NCPUS # this is the number of workers we want
                        export JULIA_NUM_THREADS=$NCPUS

                        echo running durian music after $PSEUDODEPENDS
                        sbatchid=$(sbatch \
                        --account=$SLURMACCT \
                        --partition=$SLURMPARTITION \
                        --cpus-per-task=$NCPUS \
                        --time=$slurmtimelimit \
                        --mem-per-cpu=$MEMP \
                        --job-name=$SBATCHJOBNAME \
                        --error=$SBATCHERRDIR/err_%x_%A_%a.log \
                        --out=$SBATCHOUTDIR/out_%x_%A_%a.log \
                        $SBATCHSUB | cut -f 4 -d' ')
                        NONURSMDEPENDS+=":${sbatchid}"
                        MULTISETDEPENDS+=":${sbatchid}"

                        sleep $nsleepfit

                        ######################################################################################
                        # fit durian lda
                        ######################################################################################

                        export IMPUTE_METHOD=DURIAN
                        export DECONVMETHOD=dsLDA
                        export SIMREP=42

                        SBATCHJOBNAME=${PREFIXTUPLE}_${IMPUTE_METHOD}.${DECONVMETHOD}_$suffix
                        SBATCHOUTDIR=${OUTBASEDIR}/output_logs/pseudo_fit_$IMPUTE_METHOD.$DECONVMETHOD/out
                        SBATCHERRDIR=${OUTBASEDIR}/output_logs/pseudo_fit_$IMPUTE_METHOD.$DECONVMETHOD/err
                        mkdir -p $SBATCHERRDIR
                        mkdir -p $SBATCHOUTDIR

                        export LDARUNQC=FALSE
                        export LDASCALEBLK=lognorm #lognorm
                        export LDASCALESC=column #
                        export LDASCALEFACBLK=10000
                        export MINCELLSTOPICCORP=1
                        export MCNPARTITIONS=$((NCPUS-1))
                        export MCNCHAINS=2
                        export MCTHINNING=1
                        export MCBURNRATE=0.5
                        export NJULIACORES=$((NCPUS-1)) # this should be leq the number of PBULKS 

                        SBATCHSUB=$BASEDIR/application_scripts/pseudo_array_task.sub
                        export JOBSCRIPT=$BASEDIR/application_scripts/run_imputation_methods.R
                        export nCoresAvail=$NCPUS # this is the number of workers we want
                        export JULIA_NUM_THREADS=$NCPUS

                        echo running durian dslda after $PSEUDODEPENDS
                        sbatchid=$(sbatch \
                        --account=$SLURMACCT \
                        --partition=$SLURMPARTITION \
                        --cpus-per-task=$NCPUS \
                        --time=$slurmtimelimit \
                        --mem-per-cpu=$MEMP \
                        --job-name=$SBATCHJOBNAME \
                        --error=$SBATCHERRDIR/err_%x_%A_%a.log \
                        --out=$SBATCHOUTDIR/out_%x_%A_%a.log \
                        $SBATCHSUB | cut -f 4 -d' ')
                        NONURSMDEPENDS+=":${sbatchid}"
                        MULTISETDEPENDS+=":${sbatchid}"
                done

        done
done