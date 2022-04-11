#!/bin/bash
# for he

######################################################################################
# set project environment for R and julia
######################################################################################
nsleepsim=30 # amount of time to sleep after generating sc simulation (prevent pseudo from erroring upon creation)
nsleepfit=5 # amount of time to sleep in between steps that seem to miss key environment vars
nsleeploop=5 # how long to sleep between loop iterations
dsname=He
suffix=OuterMetricsResampleShort
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
export MAXREP=20
export SUBTARGETSIZE=500
export SUBMINCELLS=10
export SUBGENERATE=0.01


MEMP=16000M # memory in mb, try increasing if nodes are not avail
export SOURCEPATH=$BASEDIR/${dsname}/durian_data
export NVME_NODESEXCLUDE=$BASEDIR/slurm_housekeeping/nodes_exclude_avx_${SLURMPARTITION}.txt

export JULIA_PROJECT=${PROJECTDIR} # make sure all workers can access the project enviroment
export JULIA_HOME=/opt/apps/julia/1.6.0/bin # make sure all workers can access the project enviroment
export R_HOME=/opt/apps/R/4.0.4/lib64/R # make sure JuliaCall/RCall can access R
export R_LIBS_USER=/data/homezvol2/mkarikom/R/x86_64-pc-linux-gnu-library/4.0 # norm() error
export LOCAL_R_LIBS_USER=/tmp/mkarikom/mylibs
export JULIA_GR_PROVIDER=GR


# export JULIA_DEPOT_PATH=/dfs5/bio/mkarikom/Julia6.0_Depot:$JULIA_DEPOT_PATH
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/dfs6/pub/mkarikom/binaries/anaconda3/lib/ # needed for pypolyagamma
export PYTHONPATH=/dfs5/bio/mkarikom/Python2.7_Pip_Packages

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

SCPARAMS=( "1,1e-6,1e-4" )
DPARAMS=( "1,1e-6,1e-4" )

export durianEps=1e-3

TypeList=( "" )
PREFIXTUPLES=( "HeLS.sense0.1.cpm.sub5000;SuarezLS.sense0.1.cpm" "HeNL.sense0.1.cpm.sub5000;SuarezNL.sense0.1.cpm" )
export PREFIXDELIM="HeLS.sense0.1.cpm.sub5000;SuarezLS.sense0.1.cpm:HeNL.sense0.1.cpm.sub5000;SuarezNL.sense0.1.cpm"

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

for SIMREP in $(seq 1 $MAXREP); do
        export SIMREP
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
                        export JOBSCRIPT=$BASEDIR/application_scripts/run_imputation_methods_clusterMetrics_outerStats_clValid_subsample.R

                        for IMPUTE_METHOD in "${IMPUTE_METHODS[@]}"; do
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
                                        --error=$SBATCHERRDIR/err_%x_%A_%a.txt \
                                        --out=$SBATCHOUTDIR/out_%x_%A_%a.txt \
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

                                SBATCHJOBNAME=${PREFIXTUPLE}_${IMPUTE_METHOD}_$suffix
                                SBATCHOUTDIR=${OUTBASEDIR}/output_logs/pseudo_fit_$IMPUTE_METHOD/out
                                SBATCHERRDIR=${OUTBASEDIR}/output_logs/pseudo_fit_$IMPUTE_METHOD/err
                                mkdir -p $SBATCHERRDIR
                                mkdir -p $SBATCHOUTDIR

                                export NJULIACORES=$((NBULK+1)) # this should be leq the number of PBULKS 
                                SBATCHSUB=$BASEDIR/application_scripts/pseudo_array_task.sub
                                export JOBSCRIPT=$BASEDIR/application_scripts/run_imputation_methods_clusterMetrics_outerStats_clValid_subsample.R
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
                                --error=$SBATCHERRDIR/err_%x_%A_%a.txt \
                                --out=$SBATCHOUTDIR/out_%x_%A_%a.txt \
                                $SBATCHSUB | cut -f 4 -d' ')
                                NONURSMDEPENDS+=":${sbatchid}"
                                MULTISETDEPENDS+=":${sbatchid}"

                                sleep $nsleepfit

                                # ######################################################################################
                                # # fit durian lda
                                # ######################################################################################

                                # export IMPUTE_METHOD=DURIAN
                                # export DECONVMETHOD=dsLDA

                                # SBATCHJOBNAME=${PREFIXTUPLE}_${IMPUTE_METHOD}_$suffix
                                # SBATCHOUTDIR=${OUTBASEDIR}/output_logs/pseudo_fit_$IMPUTE_METHOD/out
                                # SBATCHERRDIR=${OUTBASEDIR}/output_logs/pseudo_fit_$IMPUTE_METHOD/err
                                # mkdir -p $SBATCHERRDIR
                                # mkdir -p $SBATCHOUTDIR

                                # export LDARUNQC=FALSE
                                # export LDASCALEBLK=lognorm #lognorm
                                # export LDASCALESC=column #
                                # export LDASCALEFACBLK=10000
                                # export MINCELLSTOPICCORP=1
                                # export MCNPARTITIONS=$NPBULK
                                # export MCNCHAINS=2
                                # export MCTHINNING=1
                                # export MCBURNRATE=0.5
                                # export NJULIACORES=$((NBULK+1)) # this should be leq the number of PBULKS 

                                # SBATCHSUB=$BASEDIR/application_scripts/pseudo_array_task.sub
                                # export JOBSCRIPT=$BASEDIR/application_scripts/run_imputation_methods_clusterMetrics_outerStats_clValid_subsample.R
                                # export nCoresAvail=$NCPUS # this is the number of workers we want
                                # export JULIA_NUM_THREADS=$NCPUS

                                # echo running durian dslda after $PSEUDODEPENDS
                                # sbatchid=$(sbatch \
                                # --account=$SLURMACCT \
                                # --partition=$SLURMPARTITION \
                                # --cpus-per-task=$NCPUS \
                                # --time=$slurmtimelimit \
                                # --mem-per-cpu=$MEMP \
                                # --exclude=$NODESEXCLUDE \
                                # --job-name=$SBATCHJOBNAME \
                                # --error=$SBATCHERRDIR/err_%x_%A_%a.txt \
                                # --out=$SBATCHOUTDIR/out_%x_%A_%a.txt \
                                # $SBATCHSUB | cut -f 4 -d' ')
                                # NONURSMDEPENDS+=":${sbatchid}"
                                # MULTISETDEPENDS+=":${sbatchid}"
                        done

                done
        done
done