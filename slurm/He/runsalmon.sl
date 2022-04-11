#!/bin/bash

#SBATCH -A qnie_lab
#SBATCH --partition=free
#SBATCH --job-name="salmon quant"
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=40        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=7800M
#SBATCH --time=2-00:00:00          # total run time limit (HH:MM:SS)
#SBATCH --error=/dfs5/bio/mkarikom/code/DTMwork/slurm/He/logs/err/err_%x_%j.txt
#SBATCH --output=/dfs5/bio/mkarikom/code/DTMwork/slurm/He/logs/out/out_%x_%j.txt

module purge
module load salmon/1.2.1

basedir=/dfs5/bio/mkarikom/code/DTMwork/slurm/He
datadir=${basedir}/bdata_fuse

fastqfiles=(${datadir}/*)
echo ${#fastqfiles[@]}
sleep 5
indexdir=${basedir}/salmon_index
outputdir=${basedir}/salmon_output_sbatch
mkdir -p $outputdir
for i in "${fastqfiles[@]}";do
  if [[ -f "$i" ]]
  then
    echo "\n found $i\n"
    echo "running salmon quant\n\n"
    prefix="${i##/*/}"
    echo "prefix: $prefix \n"
    salmon quant -i $indexdir \
      -l A \
      -r $i \
      -o ${outputdir}/${prefix}.subset.salmon \
      -p 40 \
      --useVBOpt \
      --seqBias \
      --validateMappings

  else
      echo "\n not found $i\n"
  fi
done