#!/bin/bash

######################################################################################
# set project environment for R and julia
######################################################################################
nsleepsim=5 # amount of time to sleep after generating sc simulation (prevent pseudo from erroring upon creation)
nsleepfit=1 # amount of time to sleep in between steps that seem to miss key environment vars
nsleeploop=1 # how long to sleep between loop iterations
export nsleepdatapath=1 # how long to sleep after creating pb data path
suffix=OuterStats_clValidStabGprob

export EMDIAG=FALSE
export NPBULK=4 # the max number of pseudobulk samples within a sim
export NSIM=50
export batchFacLoc=0.1
export batchFacScale=0.2
export deFacLoc=0.5
export deFacScale=0.5
export gProbs=0.1-0.1-0.5-0.3
export deProbs=0.3-0.3-0.3-0.3 # was 0.045
export PBTRAINRATE=0.5
export SLURMPARTITION=standard
export SLURMACCT=qnie_lab
export TMPDIR=/dfs5/bio/mkarikom/sbatch_temp
export NCPUS=5
export NCPUSLOW=2
export MEMPERCPU=6000
export NODETMP=2000 # amount of /tmp space available, if this fills up then R might crash
export slurmtimelimit=0-02:30:00
export PROJECTDIR=/dfs5/bio/mkarikom/code/DURIAN
export BASEDIR=$PROJECTDIR/slurm

export RUNMASTER=$BASEDIR/durian_pseudobulk/output/pseudobulk_$suffix,gProb_$gProbs,deProb_$deProbs,bLoc_$batchFacLoc,bScale_$batchFacScale,dLoc_$deFacLoc,dScale_$deFacScale
export OUTPUTMASTER=$RUNMASTER/output.final.splatter,n_$NSIM
export URSMSCRIPT=$PROJECTDIR/URSM/scUnif_LinuxEnv.py

export JULIA_PROJECT=${PROJECTDIR} # make sure all workers can access the project enviroment
export JULIA_HOME=/opt/apps/julia/1.6.0/bin # make sure all workers can access the project enviroment
export JULIA_GR_PROVIDER=GR
export R_HOME=/opt/apps/R/4.0.4/lib64/R # make sure JuliaCall/RCall can access R
export R_LIBS_USER=/data/homezvol2/mkarikom/R/x86_64-pc-linux-gnu-library/4.0 # norm() error
export LD_LIBRARY_PATH=/opt/apps/anaconda/2020.07/lib:$LD_LIBRARY_PATH # prevent libpng16.so error when loading julia
export PYTHONPATH=/dfs6/pub/mkarikom/Python2.7_Pip_Packages

export DURIANLIB=$BASEDIR/scrabble_helper_functions/library_scrabble_clusterMetrics_clValid.R
export ETCLIB=$BASEDIR/scrabble_helper_functions/library_other_methods.R

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

######################################################################################
# ursm params
######################################################################################
# export number_of_cell_types=11 # for gupta is 11, for baron is 14
export number_of_cell_types=4 # for gupta is 11, for baron is 14
export burn_in_length=10
export gibbs_sample_number=10

######################################################################################
# pbulk params
######################################################################################
export NCELL=1600
export NGENE=800
export SIMMETHOD=splatter
export DTYPE=experiment # experiment

DROPOUTMIDS=( "0:0:0:0:4.5:4.5:4.5:4.5" "0:0:0:0:5.5:5.5:5.5:5.5" "0:0:0:0:6.5:6.5:6.5:6.5" )
export dropShapes="-1:-1:-1:-1:-1:-1:-1:-1"

export NBATCHPB=4
export NBATCHSC=4
export NPBULK=$NBATCHPB
export SPARSITY_PARAM=dmid # how the simulated dropout is named in the summarized output, eg "lambda" for downsampled real
# save all the enviroment stuff to the project dir
pip freeze > $PROJECTDIR/requirements.req
Rscript -e "sessionInfo()" >> $PROJECTDIR/session_info.req

######################################################################################
# durian/scrabble params
######################################################################################
export nEM=5 # the number of em iterations for durian and ursm
ursmEmlimit=5 # the number of em iterations for durian and ursm
export USEIRLBA=FALSE # prevent instability for benchmarks, this will take longer
export DunIterOuter=20 # 10
export DunIterInner=20 # 10
export DunSDCIters=500000
export ScrnIterOuter=20 # 20
export ScrnIterInner=20 # 20
export ScrnSDCIters=500000 # 500000
export MCNITER=5000
export DECONVGENETHRESH="-0.001"
export SCRGENETHRESH="-0.001"
export durianEps=1e-3
export SUMMARIZEDECONV=FALSE
export ERR_IN_THRESH=1e-5 # 1e-5
export ERR_OUT_THRESH=1e-7 # 1e-7
export RUNOUTERSTATS=FALSE
export RUNSTABILITY=FALSE

SCPARAMS=( "1,1e-6,1e-4" "1e-2,1e-5,1e-5" )
DPARAMS=( "1,1e-6,1e-4" "1e-2,1e-5,1e-5" )


declare -a ALLDEPENDS=''

echo starting run


for dropMids in "${DROPOUTMIDS[@]}"; do

        export dropMids

        declare -a PSEUDODEPENDS=''
        declare -a NONURSMDEPENDS=''
        declare -a URSMDEPENDS2=''
        declare -a URSMDEPENDS3=''
        declare -a URSMDEPENDS4=''

        ######################################################################################
        # pseudobulk data
        ######################################################################################
        export PBULKBASEDIR=$OUTPUTMASTER/output_splatter,simmethod_$SIMMETHOD,nsim_${NSIM},duEps_$durianEps,duEM_$nEM,urEM_$ursmEmlimit,ncell_$NCELL,ngene_$NGENE,trainrate_$PBTRAINRATE,dmid_$dropMids,suff_$suffix
        export PBULKDIR=${PBULKBASEDIR},output_fit
        export SBATCHSUB=$BASEDIR/application_scripts/pseudo_array_task.sub
        export JOBSCRIPT=$BASEDIR/durian_pseudobulk/sim_scripts/generate_splatter_k4_path_batchdrop.R

        # save the script state
        SUMMARYSCRIPT=${PBULKBASEDIR},output_scriptState
        mkdir -p ${SUMMARYSCRIPT}
        cp $ETCLIB ${SUMMARYSCRIPT}/
        cp $DURIANLIB ${SUMMARYSCRIPT}/
        cp $BASEDIR/durian_pseudobulk/sim_scripts/* ${SUMMARYSCRIPT}/

        SBATCHJOBNAME=gen_splatter_$suffix
        SBATCHOUTDIR=${PBULKBASEDIR},output_logs/pseudo_gen/out
        SBATCHERRDIR=${PBULKBASEDIR},output_logs/pseudo_gen/err
        mkdir -p $SBATCHERRDIR
        mkdir -p $SBATCHOUTDIR

        echo start splatter gen $PBULKDIR
        sbatchid=$(sbatch \
        --account=$SLURMACCT \
        --wait \
        --tmp=$NODETMP \
        --array=1-$NSIM \
        --partition=$SLURMPARTITION \
        --cpus-per-task=$NCPUSLOW \
        --time=$slurmtimelimit \
        --mem-per-cpu=$MEMPERCPU \
        --job-name=$SBATCHJOBNAME \
        --error=$SBATCHERRDIR/err_%x_%A_%a.log \
        --out=$SBATCHOUTDIR/out_%x_%A_%a.log \
        $SBATCHSUB | cut -f 4 -d' ')
        PSEUDODEPENDS+=":${sbatchid}"

        sleep $nsleepsim # sleep while splatter data is created

        ######################################################################################
        # fit non-ursm imputation methods
        ######################################################################################

        IMPUTE_METHODS=( DrImpute dropout )
        SBATCHSUB=$BASEDIR/application_scripts/pseudo_array_task.sub
        export JOBSCRIPT=$BASEDIR/application_scripts/run_imputation_methods_clusterMetrics_outerStats_clValid.R
        export nCoresAvail=$NCPUS # this is the number of workers we want
        export JULIA_NUM_THREADS=$NCPUS

        for IMPUTE_METHOD in "${IMPUTE_METHODS[@]}"; do
                export IMPUTE_METHOD
                export SIMREP=42
                echo running $IMPUTE_METHOD after $PSEUDODEPENDS
                SBATCHJOBNAME=fit_${IMPUTE_METHOD}_$suffix
                SBATCHOUTDIR=${PBULKBASEDIR},output_logs/pseudo_fit_$IMPUTE_METHOD/out
                SBATCHERRDIR=${PBULKBASEDIR},output_logs/pseudo_fit_$IMPUTE_METHOD/err
                mkdir -p $SBATCHERRDIR
                mkdir -p $SBATCHOUTDIR
                # collect the job ids `sbatchid` in an array
                sbatchid=$(sbatch \
                --account=$SLURMACCT \
                --wait \
                --tmp=$NODETMP \
                --array=1-$NSIM \
                --partition=$SLURMPARTITION \
                --cpus-per-task=$NCPUSLOW \
                --time=$slurmtimelimit \
                --mem-per-cpu=$MEMPERCPU \
                --job-name=$SBATCHJOBNAME \
                --error=$SBATCHERRDIR/err_%x_%A_%a.log \
                --out=$SBATCHOUTDIR/out_%x_%A_%a.log \
                --dependency=afterany$PSEUDODEPENDS \
                $SBATCHSUB | cut -f 4 -d' ')
                NONURSMDEPENDS+=":${sbatchid}"
                let SIMREP=SIMREP+1 
                sleep $nsleeploop
        done

        sleep $nsleepfit

        ######################################################################################
        # fit scrabble and mtscrabble
        ######################################################################################
        IMPUTE_METHODS=( SCRABBLE mtSCRABBLE )

        for IMPUTE_METHOD in "${IMPUTE_METHODS[@]}"; do
                for SCPARAM in "${SCPARAMS[@]}"; do
                        IFS="," read -r ScrabbleAlpha ScrabbleBeta ScrabbleGamma <<< "${SCPARAM}"
                        export ScrabbleAlpha
                        export ScrabbleBeta
                        export ScrabbleGamma

                        export IMPUTE_METHOD
                        export SIMREP=42
                        echo running $IMPUTE_METHOD after $PSEUDODEPENDS
                        SBATCHJOBNAME=fit_${IMPUTE_METHOD}_$suffix
                        SBATCHOUTDIR=${PBULKBASEDIR},output_logs/pseudo_fit_$IMPUTE_METHOD/out
                        SBATCHERRDIR=${PBULKBASEDIR},output_logs/pseudo_fit_$IMPUTE_METHOD/err
                        mkdir -p $SBATCHERRDIR
                        mkdir -p $SBATCHOUTDIR
                        # collect the job ids `sbatchid` in an array
                        sbatchid=$(sbatch \
                        --account=$SLURMACCT \
                        --wait \
                        --tmp=$NODETMP \
                        --array=1-$NSIM \
                        --partition=$SLURMPARTITION \
                        --cpus-per-task=$NCPUS \
                        --time=$slurmtimelimit \
                        --mem-per-cpu=$MEMPERCPU \
                        --job-name=$SBATCHJOBNAME \
                        --error=$SBATCHERRDIR/err_%x_%A_%a.log \
                        --out=$SBATCHOUTDIR/out_%x_%A_%a.log \
                        --dependency=afterany$PSEUDODEPENDS \
                        $SBATCHSUB | cut -f 4 -d' ')
                        NONURSMDEPENDS+=":${sbatchid}"
                        let SIMREP=SIMREP+1 
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

                SBATCHJOBNAME=fit_durian_${DECONVMETHOD}_$suffix
                SBATCHOUTDIR=${PBULKBASEDIR},output_logs/pseudo_fit_durian_${DECONVMETHOD}/out
                SBATCHERRDIR=${PBULKBASEDIR},output_logs/pseudo_fit_durian_${DECONVMETHOD}/err
                mkdir -p $SBATCHERRDIR
                mkdir -p $SBATCHOUTDIR

                export INITSCRABBLE=FALSE

                export NJULIACORES=$((NPBULK+1)) # this should be leq the number of PBULKS 
                SBATCHSUB=$BASEDIR/application_scripts/pseudo_array_task.sub
                export JOBSCRIPT=$BASEDIR/application_scripts/run_imputation_methods_clusterMetrics_outerStats_clValid.R
                export nCoresAvail=$NCPUS # this is the number of workers we want
                export JULIA_NUM_THREADS=$NCPUS

                echo running durian music after $PSEUDODEPENDS
                sbatchid=$(sbatch \
                --account=$SLURMACCT \
                --wait \
                --tmp=$NODETMP \
                --array=1-$NSIM \
                --partition=$SLURMPARTITION \
                --cpus-per-task=$NCPUS \
                --time=$slurmtimelimit \
                --mem-per-cpu=$MEMPERCPU \
                --job-name=$SBATCHJOBNAME \
                --error=$SBATCHERRDIR/err_%x_%A_%a.log \
                --out=$SBATCHOUTDIR/out_%x_%A_%a.log \
                --dependency=afterany$PSEUDODEPENDS \
                $SBATCHSUB | cut -f 4 -d' ')
                NONURSMDEPENDS+=":${sbatchid}"

                sleep $nsleepfit

                ######################################################################################
                # fit durian lda
                ######################################################################################

                export IMPUTE_METHOD=DURIAN
                export DECONVMETHOD=dsLDA
                export SIMREP=42

                SBATCHJOBNAME=fit_durian_${DECONVMETHOD}_$suffix
                SBATCHOUTDIR=${PBULKBASEDIR},output_logs/pseudo_fit_durian_${DECONVMETHOD}/out
                SBATCHERRDIR=${PBULKBASEDIR},output_logs/pseudo_fit_durian_${DECONVMETHOD}/err
                mkdir -p $SBATCHERRDIR
                mkdir -p $SBATCHOUTDIR

                export INITSCRABBLE=FALSE
                export LDARUNQC=FALSE
                export LDASCALEBLK=lognorm #lognorm
                export LDASCALESC=column #
                export LDASCALEFACBLK=10000
                export MINCELLSTOPICCORP=1
                export MCNPARTITIONS=$NPBULK
                export MCNCHAINS=2
                export MCTHINNING=1
                export MCBURNRATE=0.5
                export NJULIACORES=$((NPBULK+1)) # this should be leq the number of PBULKS 

                SBATCHSUB=$BASEDIR/application_scripts/pseudo_array_task.sub
                export JOBSCRIPT=$BASEDIR/application_scripts/run_imputation_methods_clusterMetrics_outerStats_clValid.R
                export nCoresAvail=$NCPUS # this is the number of workers we want
                export JULIA_NUM_THREADS=$NCPUS

                echo running durian dslda after $PSEUDODEPENDS
                sbatchid=$(sbatch \
                --account=$SLURMACCT \
                --wait \
                --tmp=$NODETMP \
                --array=1-$NSIM \
                --partition=$SLURMPARTITION \
                --cpus-per-task=$NCPUS \
                --time=$slurmtimelimit \
                --mem-per-cpu=$MEMPERCPU \
                --job-name=$SBATCHJOBNAME \
                --error=$SBATCHERRDIR/err_%x_%A_%a.log \
                --out=$SBATCHOUTDIR/out_%x_%A_%a.log \
                --dependency=afterany$PSEUDODEPENDS \
                $SBATCHSUB | cut -f 4 -d' ')
                NONURSMDEPENDS+=":${sbatchid}"
        done


        ######################################################################################
        # pre-generate all the ursm tmp files
        ######################################################################################
        SBATCHJOBNAME=pre_ursm_$suffix
        SBATCHOUTDIR=${PBULKBASEDIR},output_logs/pseudo_pre-ursm/out
        SBATCHERRDIR=${PBULKBASEDIR},output_logs/pseudo_pre-ursm/err
        mkdir -p $SBATCHERRDIR
        mkdir -p $SBATCHOUTDIR

        SBATCHSUB=$BASEDIR/application_scripts/pseudo_array_task.sub
        export JOBSCRIPT=$BASEDIR/application_scripts/prep_ursm.R

        echo running ursm prep after $PSEUDODEPENDS
        sbatchid=$(sbatch \
        --account=$SLURMACCT \
        --wait \
        --tmp=$NODETMP \
        --array=1-$NSIM \
        --partition=$SLURMPARTITION \
        --cpus-per-task=$NCPUSLOW \
        --time=$slurmtimelimit \
        --mem-per-cpu=$MEMPERCPU \
        --job-name=$SBATCHJOBNAME \
        --error=$SBATCHERRDIR/err_%x_%A_%a.log \
        --out=$SBATCHOUTDIR/out_%x_%A_%a.log \
        --dependency=afterany$PSEUDODEPENDS \
        $SBATCHSUB | cut -f 4 -d' ')
        URSMDEPENDS2+=":${sbatchid}"

        ######################################################################################
        # fit ursm
        ######################################################################################
        SBATCHJOBNAME=fit_ursm_$suffix
        SBATCHSUB=$BASEDIR/application_scripts/pseudo_array_task_ursm.sub
        SIMREP=0

        SBATCHOUTDIR=${PBULKBASEDIR},output_logs/pseudo_fit_ursm/out
        SBATCHERRDIR=${PBULKBASEDIR},output_logs/pseudo_fit_ursm/err
        mkdir -p $SBATCHERRDIR
        mkdir -p $SBATCHOUTDIR

        export EM_maxiter=$ursmEmlimit
        export output_prefix=gemout_

        echo running ursm fit after $URSMDEPENDS2
        sbatchid=$(sbatch \
        --account=$SLURMACCT \
        --tmp=$NODETMP \
        --array=1-$NSIM \
        --partition=$SLURMPARTITION \
        --cpus-per-task=$NCPUSLOW \
        --time=$slurmtimelimit \
        --mem-per-cpu=$MEMPERCPU \
        --job-name=$SBATCHJOBNAME \
        --error=$SBATCHERRDIR/err_%x_%A_%a.log \
        --out=$SBATCHOUTDIR/out_%x_%A_%a.log \
        --dependency=afterany$URSMDEPENDS2 \
        $SBATCHSUB | cut -f 4 -d' ')
        URSMDEPENDS3+=":${sbatchid}"

        ######################################################################################
        # analyze ursm
        ######################################################################################
        SBATCHJOBNAME=err_ursm_$suffix
        SBATCHOUTDIR=${PBULKBASEDIR},output_logs/pseudo_analyze-ursm/out
        SBATCHERRDIR=${PBULKBASEDIR},output_logs/pseudo_analyze-ursm/err
        mkdir -p $SBATCHERRDIR
        mkdir -p $SBATCHOUTDIR

        SBATCHSUB=$BASEDIR/application_scripts/pseudo_array_task.sub
        export JOBSCRIPT=$BASEDIR/application_scripts/analyze_ursm.R


        echo running ursm err after $URSMDEPENDS2
        sbatchid=$(sbatch \
        --account=$SLURMACCT \
        --wait \
        --tmp=$NODETMP \
        --array=1-$NSIM \
        --partition=$SLURMPARTITION \
        --cpus-per-task=$NCPUSLOW \
        --time=$slurmtimelimit \
        --mem-per-cpu=$MEMPERCPU \
        --job-name=$SBATCHJOBNAME \
        --error=$SBATCHERRDIR/err_%x_%A_%a.log \
        --out=$SBATCHOUTDIR/out_%x_%A_%a.log \
        --dependency=afterany$URSMDEPENDS3 \
        $SBATCHSUB | cut -f 4 -d' ')
        URSMDEPENDS4+=":${sbatchid}"

        export SUMMARYDEPENDS=$NONURSMDEPENDS$URSMDEPENDS4
        # sleep $nsleepfit

        ######################################################################################
        # summarize final output of all methods
        ######################################################################################
        SBATCHJOBNAME=summarize_all_$suffix
        SBATCHSUB=$BASEDIR/application_scripts/pseudo_array_task.sub
        export JOBSCRIPT=$BASEDIR/application_scripts/summarize_plot_final_pseudo.R
        export SUMMARYFINAL=${PBULKBASEDIR},output_summaryFinal

        SBATCHOUTDIR=${PBULKBASEDIR},output_logs/pseudo_summarize_final/out
        SBATCHERRDIR=${PBULKBASEDIR},output_logs/pseudo_summarize_final/err
        mkdir -p $SBATCHERRDIR
        mkdir -p $SBATCHOUTDIR

        sbatchid=$(sbatch \
        --account=$SLURMACCT \
        --tmp=$NODETMP \
        --partition=$SLURMPARTITION \
        --cpus-per-task=$NCPUSLOW \
        --time=$slurmtimelimit \
        --mem-per-cpu=$MEMPERCPU \
        --job-name=$SBATCHJOBNAME \
        --error=$SBATCHERRDIR/err_%x_%A_%a.log \
        --out=$SBATCHOUTDIR/out_%x_%A_%a.log \
        --dependency=afterany$SUMMARYDEPENDS \
        $SBATCHSUB | cut -f 4 -d' ')
        export ALLDEPENDS=:${sbatchid}${ALLDEPENDS}

        ######################################################################################
        # summarize iter output of durian
        ######################################################################################
        SBATCHJOBNAME=summarize_iter_$suffix
        SBATCHSUB=$BASEDIR/application_scripts/pseudo_array_task.sub
        export JOBSCRIPT=$BASEDIR/application_scripts/summarize_plot_iter_pseudo.R
        export SUMMARYITER=${PBULKBASEDIR},output_summaryIter

        SBATCHOUTDIR=${PBULKBASEDIR},output_logs/pseudo_summarize_iter/out
        SBATCHERRDIR=${PBULKBASEDIR},output_logs/pseudo_summarize_iter/err
        mkdir -p $SBATCHERRDIR
        mkdir -p $SBATCHOUTDIR

        sbatchid=$(sbatch \
        --account=$SLURMACCT \
        --tmp=$NODETMP \
        --wait \
        --partition=$SLURMPARTITION \
        --cpus-per-task=$NCPUSLOW \
        --time=$slurmtimelimit \
        --mem-per-cpu=$MEMPERCPU \
        --job-name=$SBATCHJOBNAME \
        --error=$SBATCHERRDIR/err_%x_%A_%a.log \
        --out=$SBATCHOUTDIR/out_%x_%A_%a.log \
        --dependency=afterany$SUMMARYDEPENDS \
        $SBATCHSUB | cut -f 4 -d' ')
        export ALLDEPENDS=:${sbatchid}${ALLDEPENDS}
done
######################################################################################
# combine all output for fig2
######################################################################################
SBATCHJOBNAME=combine_$suffix
SBATCHSUB=$BASEDIR/application_scripts/pseudo_array_task.sub
export JOBSCRIPT=$BASEDIR/application_scripts/summarize_plot_final_pseudo_combine.R
export SUMMARYCOMBINE=${OUTPUTMASTER},output_summaryCombine

SBATCHOUTDIR=${OUTPUTMASTER},output_logs/pseudo_summarize_combine/out
SBATCHERRDIR=${OUTPUTMASTER},output_logs/pseudo_summarize_combine/err
mkdir -p $SBATCHERRDIR
mkdir -p $SBATCHOUTDIR

sbatchid=$(sbatch \
--account=$SLURMACCT \
--tmp=$NODETMP \
--partition=$SLURMPARTITION \
--cpus-per-task=$NCPUSLOW \
--time=$slurmtimelimit \
--mem-per-cpu=$MEMPERCPU \
--job-name=$SBATCHJOBNAME \
--error=$SBATCHERRDIR/err_%x_%A_%a.log \
--out=$SBATCHOUTDIR/out_%x_%A_%a.log \
--dependency=afterany$ALLDEPENDS \
$SBATCHSUB | cut -f 4 -d' ')

######################################################################################
# combine all output for fig2
######################################################################################
SBATCHJOBNAME=combine_$suffix
SBATCHSUB=$BASEDIR/application_scripts/pseudo_array_task.sub
export JOBSCRIPT=$BASEDIR/application_scripts/summarize_plot_final_pseudo_combine_param.R
export SUMMARYCOMBINE=${OUTPUTMASTER},output_summaryCombineParam

SBATCHOUTDIR=${OUTPUTMASTER},output_logs/pseudo_summarize_combine_param/out
SBATCHERRDIR=${OUTPUTMASTER},output_logs/pseudo_summarize_combine_param/err
mkdir -p $SBATCHERRDIR
mkdir -p $SBATCHOUTDIR

sbatchid=$(sbatch \
--account=$SLURMACCT \
--tmp=$NODETMP \
--partition=$SLURMPARTITION \
--cpus-per-task=$NCPUSLOW \
--time=$slurmtimelimit \
--mem-per-cpu=$MEMPERCPU \
--job-name=$SBATCHJOBNAME \
--error=$SBATCHERRDIR/err_%x_%A_%a.log \
--out=$SBATCHOUTDIR/out_%x_%A_%a.log \
--dependency=afterany$ALLDEPENDS \
$SBATCHSUB | cut -f 4 -d' ')

######################################################################################
# combine all iter output for fig2
######################################################################################
SBATCHJOBNAME=iter_combine_$suffix
SBATCHSUB=$BASEDIR/application_scripts/pseudo_array_task.sub
export JOBSCRIPT=$BASEDIR/application_scripts/summarize_plot_iter_pseudo_combine.R
export SUMMARYCOMBINE_ITER=${OUTPUTMASTER},output_summaryCombineIter

SBATCHOUTDIR=${OUTPUTMASTER},output_logs/pseudo_summarize_combine_iter/out
SBATCHERRDIR=${OUTPUTMASTER},output_logs/pseudo_summarize_combine_iter/err
mkdir -p $SBATCHERRDIR
mkdir -p $SBATCHOUTDIR

sbatchid=$(sbatch \
--account=$SLURMACCT \
--tmp=$NODETMP \
--partition=$SLURMPARTITION \
--cpus-per-task=$NCPUSLOW \
--time=$slurmtimelimit \
--mem-per-cpu=$MEMPERCPU \
--job-name=$SBATCHJOBNAME \
--error=$SBATCHERRDIR/err_%x_%A_%a.log \
--out=$SBATCHOUTDIR/out_%x_%A_%a.log \
--dependency=afterany$ALLDEPENDS \
$SBATCHSUB | cut -f 4 -d' ')