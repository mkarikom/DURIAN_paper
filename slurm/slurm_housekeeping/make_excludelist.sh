#!/bin/bash

ls $NVME_NODEDIR > $NVME_NODEFILE
sleep 5

ls -lavh $ALLNODES
ls -lavh $NVME_NODEFILE

# make the exclude list
awk 'NR==FNR {key[$1]; next} !($1 in key)' $NVME_NODEFILE $ALLNODES
awk 'NR==FNR {key[$1]; next} !($1 in key)' $NVME_NODEFILE $ALLNODES > $NVME_NODESEXCLUDE

ls -alvh $NVME_NODESEXCLUDE