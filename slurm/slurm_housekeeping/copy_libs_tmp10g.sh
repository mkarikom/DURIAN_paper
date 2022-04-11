#!/bin/bash

######################################################################################
# copy r libs to file and make exclude list
######################################################################################

nodeprefs=( "hpc3-14-" "hpc3-15-" "hpc3-17-" "hpc3-20-" "hpc3-21-" "hpc3-22-" )
nodeinds=({00..48})
suffix=housekeeping
NCPUS=1

rm -fr /dfs5/bio/mkarikom/temp/DURIAN/slurm/slurm_housekeeping/output_logs_tmp10g
rm -fr /dfs5/bio/mkarikom/temp/DURIAN/slurm/slurm_housekeeping/nodes_ready_tmp10g
rm /dfs5/bio/mkarikom/temp/DURIAN/slurm/slurm_housekeeping/nodes_ready_tmp10g.txt
# rm /dfs5/bio/mkarikom/temp/DURIAN/slurm/slurm_housekeeping/nodes_exclude.txt

export SLURMPARTITION=standard
export SLURMTIMELIMIT=0-01:00:00
export PROJECTDIR=/dfs6/pub/mkarikom/code/DURIAN_paper_clean
export BASEDIR=$PROJECTDIR/slurm
export NVME_NODEFILE=$BASEDIR/slurm_housekeeping/nodes_ready_tmp10g.txt
export NVME_NODEDIR=$BASEDIR/slurm_housekeeping/nodes_ready_tmp10g
export NVME_NODESEXCLUDE=$BASEDIR/slurm_housekeeping/nodes_exclude_tmp10g.txt
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
SBATCHJOBNAME=tmp10g_clb
SBATCHOUTDIR=$BASEDIR/slurm_housekeeping/output_logs_tmp10g/out
SBATCHERRDIR=$BASEDIR/slurm_housekeeping/output_logs_tmp10g/err
mkdir -p $SBATCHERRDIR
mkdir -p $SBATCHOUTDIR

# save all node names
export ALLNODES=$BASEDIR/slurm_housekeeping/nodes_all_tmp10g.txt
sinfo -o "%n"| awk '{if (NR!=1) {print $1}}'|less > $ALLNODES

mkdir -p $NVME_NODEDIR

declare -a NODEJOBS=''

for nodeprefix in "${nodeprefs[@]}"; do
    for ind in "${nodeinds[@]}"; do
        echo copying libs to node ${nodeprefix}${ind}
        sbatchid=$(sbatch \
        --account=qnie_lab \
        --partition=$SLURMPARTITION \
        --cpus-per-task=$NCPUS \
        --tmp=10000M \
        --nodelist=${nodeprefix}${ind} \
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
done

echo nodejobs $NODEJOBS

SBATCHSUB=$BASEDIR/slurm_housekeeping/make_excludelist.sub
SBATCHJOBNAME=tmp10g_nodelist
sbatch \
--account=qnie_lab \
--partition=$SLURMPARTITION \
--cpus-per-task=$NCPUS \
--constraint=nvme \
--time=$SLURMTIMELIMIT \
--mem-per-cpu=$MEMPERCPU \
--job-name=$SBATCHJOBNAME \
--dependency=afterany$NODEJOBS \
--error=$SBATCHERRDIR/err_%x_%A_%a.log \
--out=$SBATCHOUTDIR/out_%x_%A_%a.log \
$SBATCHSUB