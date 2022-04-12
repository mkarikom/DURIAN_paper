#!/bin/bash
# for baron

######################################################################################
# set project environment for R and julia
######################################################################################
nsleepsim=30 # amount of time to sleep after generating sc simulation (prevent pseudo from erroring upon creation)
nsleepfit=5 # amount of time to sleep in between steps that seem to miss key environment vars
nsleeploop=5 # how long to sleep between loop iterations
dsname=Baron
suffix=OuterMetrics1K05
prefix=$dsname
SLURMACCT=qnie_lab

export EMDIAG=FALSE
export SLURMPARTITION=standard
export slurmtimelimit=3-00:00:00
export PROJECTDIR=/dfs6/pub/mkarikom/code/DURIAN_paper_clean
export BASEDIR=$PROJECTDIR/slurm
export OUTPUTMASTER=$BASEDIR/${dsname}/output.clusterMetrics.$SLURMPARTITION.$suffix
export NBULK=7
export NCPUS=10

MEMP=16000M # memory in mb, try increasing if nodes are not avail
export SOURCEPATH=$BASEDIR/${dsname}/durian_data
export NVME_NODESEXCLUDE=$BASEDIR/slurm_housekeeping/nodes_exclude_avx_${SLURMPARTITION}.txt

export JULIA_PROJECT=${PROJECTDIR} # make sure all workers can access the project enviroment
export JULIA_HOME=/opt/apps/julia/1.6.0/bin # make sure all workers can access the project enviroment
export R_HOME=/opt/apps/R/4.0.4/lib64/R # make sure JuliaCall/RCall can access R
export R_LIBS_USER=/data/homezvol2/mkarikom/R/x86_64-pc-linux-gnu-library/4.0 # norm() error
export JULIA_GR_PROVIDER=GR
export LD_LIBRARY_PATH=/opt/apps/anaconda/2020.07/lib:$LD_LIBRARY_PATH # prevent libpng16.so error when loading julia
export PYTHONPATH=/dfs6/pub/mkarikom/Python2.7_Pip_Packages

module purge
module load zlib
# module load eigen
# module load hdf5
# module load gcc
# module load lapack
# module load OpenBLAS
# module load mkl
module load julia/1.6.0
module load R/4.0.4
module load python/2.7.17 # needed for ursm, pypolyagamma

export DURIANLIB=$BASEDIR/scrabble_helper_functions/library_scrabble_clusterMetrics_clValid.R
export ETCLIB=$BASEDIR/scrabble_helper_functions/library_other_methods.R
export SIGNALINGLIB=$BASEDIR/scrabble_helper_functions/library_signaling.R
export MULTISETLIB=$BASEDIR/scrabble_helper_functions/library_signaling_multiset.R

# save all the enviroment stuff to the project dir
pip freeze > $PROJECTDIR/requirements.txt
Rscript -e "sessionInfo()" >> $PROJECTDIR/session_info.txt

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

export ERR_IN_THRESH=1e-5 # 1e-5
export ERR_OUT_THRESH=1e-7 # 1e-7
export RUNOUTERSTATS=TRUE
export RUNSTABILITY=FALSE
export INITSCRABBLE=FALSE
export SUMMARIZEDECONV=TRUE

export DECONVGENETHRESH=-0.001
export SCRGENETHRESH=-0.001

SCPARAMS=( "1,1e-6,1e-4" "1e-2,1e-5,1e-5" )
DPARAMS=( "1,1e-6,1e-4" "1e-2,1e-5,1e-5" )

export durianEps=1e-3

TypeList=( "" )
PREFIXTUPLES=( "BaronSC.DM.isletVST1K05;SegerstolpeBulk.DM.cpm" "BaronSC.H.isletVST1K05;SegerstolpeBulk.H.cpm" )
export PREFIXDELIM="BaronSC.DM.isletVST1K05;SegerstolpeBulk.DM.cpm:BaronSC.H.isletVST1K05;SegerstolpeBulk.H.cpm"

export PAIREDPRIORITY="BaronSC.H"
export COMMONSUFF=".isletVST1K05"

export CELLCHATDB="CellChatDB.human"
nsleepsim=60 # amount of time to sleep after generating sc simulation (prevent pseudo from erroring upon creation)
nsleepfit=2 # amount of time to sleep in between steps that seem to miss key environment vars
nsleeploop=0.05 # how long to sleep between loop iterations
export nsleepdatapath=1 # how long to sleep after creating pb data path
#################################################################


######################################################################################
# ursm params
######################################################################################
# export number_of_cell_types=11 # for gupta is 11, for baron is 14
export number_of_cell_types=14 # for gupta is 11, for baron is 14
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
                IMPUTE_METHODS=( DrImpute dropout )
                export JOBSCRIPT=$BASEDIR/application_scripts/run_imputation_methods_clusterMetrics_outerStats_clValid.R

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
                        --exclude=$NVME_NODESEXCLUDE \
                        --job-name=$SBATCHJOBNAME \
                        --error=$SBATCHERRDIR/err_%x_%A.txt \
                        --out=$SBATCHOUTDIR/out_%x_%A.txt \
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
                                --exclude=$NODESEXCLUDE \
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

                        SBATCHJOBNAME=${PREFIXTUPLE}_${IMPUTE_METHOD}_$suffix
                        SBATCHOUTDIR=${OUTBASEDIR}/output_logs/pseudo_fit_$IMPUTE_METHOD/out
                        SBATCHERRDIR=${OUTBASEDIR}/output_logs/pseudo_fit_$IMPUTE_METHOD/err
                        mkdir -p $SBATCHERRDIR
                        mkdir -p $SBATCHOUTDIR

                        export NJULIACORES=$((NBULK+1)) # this should be leq the number of PBULKS 
                        SBATCHSUB=$BASEDIR/application_scripts/pseudo_array_task.sub
                        export JOBSCRIPT=$BASEDIR/application_scripts/run_imputation_methods_clusterMetrics_outerStats_clValid.R
                        export nCoresAvail=$NCPUS # this is the number of workers we want
                        export JULIA_NUM_THREADS=$NCPUS

                        echo running durian music after $PSEUDODEPENDS
                        sbatchid=$(sbatch \
                        --account=$SLURMACCT \
                        --partition=$SLURMPARTITION \
                        --cpus-per-task=$NCPUS \
                        --time=$slurmtimelimit \
                        --mem-per-cpu=$MEMP \
                        --exclude=$NODESEXCLUDE \
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

                        SBATCHJOBNAME=${PREFIXTUPLE}_${IMPUTE_METHOD}_$suffix
                        SBATCHOUTDIR=${OUTBASEDIR}/output_logs/pseudo_fit_$IMPUTE_METHOD/out
                        SBATCHERRDIR=${OUTBASEDIR}/output_logs/pseudo_fit_$IMPUTE_METHOD/err
                        mkdir -p $SBATCHERRDIR
                        mkdir -p $SBATCHOUTDIR

                        export LDARUNQC=FALSE
                        export LDASCALEBLK=lognorm #lognorm
                        export LDASCALESC=column #
                        export LDASCALEFACBLK=10000
                        export MINCELLSTOPICCORP=1
                        export MCNPARTITIONS=$NPBULK
                        export MCNCHAINS=2
                        export MCTHINNING=1
                        export MCBURNRATE=0.5
                        export NJULIACORES=$((NBULK+1)) # this should be leq the number of PBULKS 

                        SBATCHSUB=$BASEDIR/application_scripts/pseudo_array_task.sub
                        export JOBSCRIPT=$BASEDIR/application_scripts/run_imputation_methods_clusterMetrics_outerStats_clValid.R
                        export nCoresAvail=$NCPUS # this is the number of workers we want
                        export JULIA_NUM_THREADS=$NCPUS

                        echo running durian dslda after $PSEUDODEPENDS
                        sbatchid=$(sbatch \
                        --account=$SLURMACCT \
                        --partition=$SLURMPARTITION \
                        --cpus-per-task=$NCPUS \
                        --time=$slurmtimelimit \
                        --mem-per-cpu=$MEMP \
                        --exclude=$NODESEXCLUDE \
                        --job-name=$SBATCHJOBNAME \
                        --error=$SBATCHERRDIR/err_%x_%A_%a.log \
                        --out=$SBATCHOUTDIR/out_%x_%A_%a.log \
                        $SBATCHSUB | cut -f 4 -d' ')
                        NONURSMDEPENDS+=":${sbatchid}"
                        MULTISETDEPENDS+=":${sbatchid}"
                done

        done

        SBATCHJOBNAME=clust
        SBATCHOUTDIR=${OUTPUTMASTER}/output_logs/multisetClust/out
        SBATCHERRDIR=${OUTPUTMASTER}/output_logs/multisetClust/err
        mkdir -p $SBATCHERRDIR
        mkdir -p $SBATCHOUTDIR

        export JOBSCRIPT=$BASEDIR/application_scripts/run_clusterMetrics_real.R
        export SAVEPATH=${OUTPUTMASTER}/clusterMetrics_multiset

        echo running multiset clustermetrics
        sbatchid=$(sbatch \
        --account=$SLURMACCT \
        --partition=$SLURMPARTITION \
        --cpus-per-task=$NCPUS \
        --time=$slurmtimelimit \
        --mem-per-cpu=$MEMP \
        --exclude=$NVME_NODESEXCLUDE \
        --job-name=$SBATCHJOBNAME \
        --dependency=afterany$MULTISETDEPENDS \
        --error=$SBATCHERRDIR/err_%x_%A.txt \
        --out=$SBATCHOUTDIR/out_%x_%A.txt \
        $SBATCHSUB | cut -f 4 -d' ')
done