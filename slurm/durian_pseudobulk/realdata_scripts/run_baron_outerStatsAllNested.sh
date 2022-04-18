#!/bin/bash

######################################################################################
# set project environment for R and julia
######################################################################################
nsleepsim=120 # amount of time to sleep after generating sc simulation (prevent pseudo from erroring upon creation)
nsleepfit=5 # amount of time to sleep in between steps that seem to miss key environment vars
nsleeploop=5 # how long to sleep between loop iterations
export nsleepdatapath=1 # how long to sleep after creating pb data path
suffix=BaronOuterStatsAllNested

export EMDIAG=FALSE
export NPBULK=3 # the max number of pseudobulk samples within a sim
# export NSIM=50
export NSIM=2
export PBTRAINRATE=0.5
# export SLURMPARTITION=standard
export SLURMPARTITION=free
export SLURMACCT=qnie_lab
MEMPERCPU=10000
# export slurmtimelimit=0-02:30:00
export slurmtimelimit=0-00:30:00
export PROJECTDIR=/dfs6/pub/mkarikom/code/DURIAN_paper_clean
export BASEDIR=$PROJECTDIR/slurm
export NCPUS=$((NPBULK+1))
export NCPUSLOW=2

export RUNMASTER=$BASEDIR/durian_pseudobulk_$suffix
export URSMSCRIPT=$PROJECTDIR/URSM/scUnif_LinuxEnv.py

export JULIA_PROJECT=${PROJECTDIR} # make sure all workers can access the project enviroment
export JULIA_HOME=/opt/apps/julia/1.6.0/bin # make sure all workers can access the project enviroment
export JULIA_GR_PROVIDER=GR
export R_HOME=/opt/apps/R/4.0.4/lib64/R # make sure JuliaCall/RCall can access R
export R_LIBS_USER=/data/homezvol2/mkarikom/R/x86_64-pc-linux-gnu-library/4.0 # norm() error
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

######################################################################################
# ursm params
######################################################################################
# export number_of_cell_types=11 # for gupta is 11, for baron is 14
export number_of_cell_types=5 # for gupta is 11, for baron is 14
export burn_in_length=10
export gibbs_sample_number=10
######################################################################################
# pbulk params
######################################################################################
export NZTHRESH=0.01
export NCELL=1600 # number of random cells to select (seed is PSEUDOCOUNT), this is before train/test split
export NGENE=1000 # number of random genes to select (seed is PSEUDOCOUNT)
LAMBDAS=( 1e-5 5e-6 1e-6 )
TypeList=( alpha-beta-delta-epsilon-gamma )
export SCDATADIR=$PROJECTDIR/slurm/Baron/durian_data
export SCPREFIX="BaronSC.H_"
export SPARSITY_PARAM=lambda # how the simulated dropout is named in the summarized output, eg "lambda" for downsampled real

# save all the enviroment stuff to the project dir
pip freeze > $PROJECTDIR/requirements.req
Rscript -e "sessionInfo()" >> $PROJECTDIR/session_info.req
######################################################################################
# durian/scrabble params
######################################################################################
export nEM=5 # the number of em iterations for durian and ursm
export ursmEmlimit=5 # the number of em iterations for durian and ursm
export USEIRLBA=FALSE # prevent instability for benchmarks, this will take longer
export DunIterOuter=20
export DunIterInner=20
export DunSDCIters=500000
export ScrnIterOuter=20
export ScrnIterInner=20
export ScrnSDCIters=500000
export MCNITER=5000
export DECONVGENETHRESH=0.005
export SCRGENETHRESH="-0.001"
export durianEps=1e-3
export SUMMARIZEDECONV=FALSE
export ERR_IN_THRESH=1e-5 # 1e-5
export ERR_OUT_THRESH=1e-7 # 1e-7
export RUNOUTERSTATS=FALSE
export RUNSTABILITY=FALSE

SCPARAMS=( "1,1e-6,1e-4" "1e-2,1e-5,1e-5" )
DPARAMS=( "1,1e-6,1e-4" "1e-2,1e-5,1e-5" )


#################################################################

declare -a ALLDEPENDS=''

for SUBSETCELLTYPES in "${TypeList[@]}"; do
        for LAMBDA in "${LAMBDAS[@]}"; do
                IFS="," read -r DurianAlpha DurianBeta DurianGamma <<< "${DPARAM}"
                export DurianAlpha
                export DurianBeta
                export DurianGamma
                export LAMBDA
                export OUTPUTMASTER=$BASEDIR/durian_pseudobulk/output/durian_pseudobulk_$suffix/output.final.baron,n_$NSIM
                export SUBSETCELLTYPES

                declare -a PSEUDODEPENDS=''
                declare -a NONURSMDEPENDS=''
                declare -a URSMDEPENDS2=''
                declare -a URSMDEPENDS3=''
                declare -a URSMDEPENDS4=''

                ######################################################################################
                # pseudobulk data
                ######################################################################################
                export PBULKBASEDIR=$OUTPUTMASTER/output_pseudo,nsim_${NSIM},duEps_$durianEps,duEM_$nEM,urEM_$ursmEmlimit,nz_$NZTHRESH,ncell_$NCELL,ngene_$NGENE,trainrate_$PBTRAINRATE,lambda_$LAMBDA,suff_$suffix
                export PBULKDIR=${PBULKBASEDIR},output_fit
                export SBATCHSUB=$BASEDIR/application_scripts/pseudo_array_task.sub
                export JOBSCRIPT=$BASEDIR/durian_pseudobulk/realdata_scripts/make_pseudo_donordrop.R

                SBATCHJOBNAME=gen_pseudo_$suffix
                SBATCHOUTDIR=${PBULKBASEDIR},output_logs/pseudo_gen/out
                SBATCHERRDIR=${PBULKBASEDIR},output_logs/pseudo_gen/err
                mkdir -p $SBATCHERRDIR
                mkdir -p $SBATCHOUTDIR

                echo sched pseudo
                sbatchid=$(sbatch \
                --account=$SLURMACCT \
                --array=1-$NSIM \
                --partition=$SLURMPARTITION \
                --cpus-per-task=$NCPUS \
                --time=$slurmtimelimit \
                --mem-per-cpu=$MEMPERCPU \
                --job-name=$SBATCHJOBNAME \
                --error=$SBATCHERRDIR/err_%x_%A_%a.log \
                --out=$SBATCHOUTDIR/out_%x_%A_%a.log \
                --wait \
                $SBATCHSUB | cut -f 4 -d' ')
                PSEUDODEPENDS+=":${sbatchid}"

                sleep $nsleepsim # wait for the symsim directories to be created

                ######################################################################################
                # fit non-ursm imputation methods
                ######################################################################################


                IMPUTE_METHODS=( DrImpute dropout ALRA G2S3 CMFImpute )
                SBATCHSUB=$BASEDIR/application_scripts/pseudo_array_task.sub
                export JOBSCRIPT=$BASEDIR/application_scripts/run_imputation_methods.R
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
                        export JOBSCRIPT=$BASEDIR/application_scripts/run_imputation_methods.R
                        export nCoresAvail=$NCPUS # this is the number of workers we want
                        export JULIA_NUM_THREADS=$NCPUS

                        echo running durian music after $PSEUDODEPENDS
                        sbatchid=$(sbatch \
                        --account=$SLURMACCT \
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
                        export JOBSCRIPT=$BASEDIR/application_scripts/run_imputation_methods.R
                        export nCoresAvail=$NCPUS # this is the number of workers we want
                        export JULIA_NUM_THREADS=$NCPUS

                        echo running durian dslda after $PSEUDODEPENDS
                        sbatchid=$(sbatch \
                        --account=$SLURMACCT \
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
        --partition=$SLURMPARTITION \
        --cpus-per-task=$NCPUSLOW \
        --time=$slurmtimelimit \
        --mem-per-cpu=$MEMPERCPU \
        --job-name=$SBATCHJOBNAME \
        --error=$SBATCHERRDIR/err_%x_%A_%a.log \
        --out=$SBATCHOUTDIR/out_%x_%A_%a.log \
        --dependency=afterany$ALLDEPENDS \
        $SBATCHSUB | cut -f 4 -d' ')

done