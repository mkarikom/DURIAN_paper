#!/bin/bash

######################################################################################
# copy r libs to file and make exclude list
# to cancel stalled run `scancel -u mkarikom -A qnie_lab -p free`
######################################################################################

sinfo -N | awk '{if (NR!=1) {print $1} }' > ./slurm/slurm_housekeeping/temp.txt
readarray -t nodeprefs0 < ./slurm/slurm_housekeeping/temp.txt
nodeprefs=($(echo "${nodeprefs0[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' '))

suffix=housekeeping
NCPUS=1

rm -fr /dfs5/bio/mkarikom/temp/DURIAN/slurm/slurm_housekeeping/output_logs_avx
rm -fr /dfs5/bio/mkarikom/temp/DURIAN/slurm/slurm_housekeeping/nodes_ready_avx
rm /dfs5/bio/mkarikom/temp/DURIAN/slurm/slurm_housekeeping/nodes_ready_avx.txt
# rm /dfs5/bio/mkarikom/temp/DURIAN/slurm/slurm_housekeeping/nodes_exclude.txt

export SLURMPARTITION=free
export SLURMTIMELIMIT=0-00:10:00
export PROJECTDIR=/dfs6/pub/mkarikom/code/DURIAN_paper_clean
export BASEDIR=$PROJECTDIR/slurm
export NVME_NODEFILE=$BASEDIR/slurm_housekeeping/nodes_ready_avx.txt
export NVME_NODEDIR=$BASEDIR/slurm_housekeeping/nodes_ready_avx
export NVME_NODESEXCLUDE=$BASEDIR/slurm_housekeeping/nodes_exclude_avx_$SLURMPARTITION.txt
export MEMPERCPU=500M

export R_LIBS_USER=/data/homezvol2/mkarikom/R/x86_64-pc-linux-gnu-library/4.0 # norm() error
# export R_LIBS_USER=/dfs5/bio/mkarikom/R4.0_Packages # irlba() error
export JULIA_GR_PROVIDER=GR

# `export LOCAL_R_LIBS_USER=NULL` to use /data/homezvol2/mkarikom/R/x86_64-pc-linux-gnu-library/4.0
export LOCAL_R_LIBS_USER=/tmp/mkarikom/mylibs
export SOURCE_LIB_ZIP=mylibs.tar.gz
module purge
module load zlib
module load R/4.0.4

SBATCHSUB=$BASEDIR/slurm_housekeeping/transfer.sub
SBATCHJOBNAME=avx_clb
SBATCHOUTDIR=$BASEDIR/slurm_housekeeping/output_logs_avx/out
SBATCHERRDIR=$BASEDIR/slurm_housekeeping/output_logs_avx/err
mkdir -p $SBATCHERRDIR
mkdir -p $SBATCHOUTDIR

# save all node names
export ALLNODES=$BASEDIR/slurm_housekeeping/nodes_all_avx.txt
sinfo -o "%n"| awk '{if (NR!=1) {print $1}}'|less > $ALLNODES

mkdir -p $NVME_NODEDIR

declare -a NODEJOBS=''

for nodeprefix in "${nodeprefs[@]}"; do
    echo copying libs to node ${nodeprefix}
    sbatchid=$(sbatch \
    --account=qnie_lab \
    --partition=$SLURMPARTITION \
    --cpus-per-task=$NCPUS \
    --nodelist=${nodeprefix} \
    --time=$SLURMTIMELIMIT \
    --mem-per-cpu=$MEMPERCPU \
    --job-name=$SBATCHJOBNAME \
    --error=$SBATCHERRDIR/err_%x_%A_%a.log \
    --out=$SBATCHOUTDIR/out_%x_%A_%a.log \
    $SBATCHSUB | cut -f 4 -d' ')
    echo sbatch id is $sbatchid
    if [ ${#sbatchid} -ge 2 ];
    then
        NODEJOBS+=":${sbatchid}"
    fi
done

echo nodejobs $NODEJOBS

SBATCHSUB=$BASEDIR/slurm_housekeeping/make_excludelist.sub
SBATCHJOBNAME=avx_nodelist
sbatch \
--account=qnie_lab \
--partition=debug \
--cpus-per-task=$NCPUS \
--time=$SLURMTIMELIMIT \
--mem-per-cpu=$MEMPERCPU \
--job-name=$SBATCHJOBNAME \
--dependency=afterany$NODEJOBS \
--error=$SBATCHERRDIR/err_%x_%A_%a.log \
--out=$SBATCHOUTDIR/out_%x_%A_%a.log \
$SBATCHSUB