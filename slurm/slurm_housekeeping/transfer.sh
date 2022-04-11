#!/bin/bash

####################################################
# select r libs user
####################################################

# assumes we already ran `tar -zcvf ~/rlibs/mylibs.tar.gz ~/R/x86_64-pc-linux-gnu-library/4.0`
# check if caching is desired
echo job id $SLURM_JOB_ID
echo on node $SLURM_JOB_NODELIST
rm -fr /tmp/mkarikom/TRUE
rm -fr /tmp/mkarikom/*
sleep 5
if [ ~/rlibs/$SOURCE_LIB_ZIP -nt /tmp/mkarikom/$SOURCE_LIB_ZIP ] || [ ! -f /tmp/mkarikom/$SOURCE_LIB_ZIP ]
then
    mkdir -p /tmp/mkarikom
    echo "removing outdated /tmp/mkarikom/$SOURCE_LIB_ZIP"
    rm /tmp/mkarikom/$SOURCE_LIB_ZIP
    echo "removing outdated $LOCAL_R_LIBS_USER"
    rm -fr $LOCAL_R_LIBS_USER
    echo "transferring ~/rlibs/$SOURCE_LIB_ZIP"
    rsync -Prv ~/rlibs/$SOURCE_LIB_ZIP /tmp/mkarikom/
    cd /tmp/mkarikom
    echo "unpacking R libs"
    tar xzf $SOURCE_LIB_ZIP
    echo "renaming local libs directory"
    mv /tmp/mkarikom/data/homezvol2/mkarikom/R/x86_64-pc-linux-gnu-library/4.0 $LOCAL_R_LIBS_USER
    rm -fr /tmp/mkarikom/data
    export R_LIBS_USER=$LOCAL_R_LIBS_USER
    echo using $R_LIBS_USER
    echo adding $NVME_NODEDIR/$SLURM_JOB_NODELIST to ready nodes
    touch $NVME_NODEDIR/$SLURM_JOB_NODELIST
else
    echo "cached libs already exist"
    export R_LIBS_USER=$LOCAL_R_LIBS_USER
    echo using $R_LIBS_USER
    echo adding $NVME_NODEDIR/$SLURM_JOB_NODELIST to ready nodes
    touch $NVME_NODEDIR/$SLURM_JOB_NODELIST
fi