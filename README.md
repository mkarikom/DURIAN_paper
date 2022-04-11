# DURIAN
This repository generates benchmarks and signaling analyses from the paper.

In order for R to call the dsLDA sampler (written in julia and available as DistributedTopicModels.jl), we need to eliminate any potential conflicts between on hpc, R (compiled by hpc admins) and julia (downloaded as a binary).

### how to install RCall
1. `rm ~/.julia` # this is destructive but prevents cryptic errors down the line
2. `rm Project.toml && rm Manifest.toml`
3. `export JULIA_PKG_SERVER=pkg.julialang.org` # this may provide a speedup as the ~/.julia environment is rebuilt
4. `export R_HOME=/opt/apps/R/4.0.4/lib64/R` # (or `module load R`) make sure rcall.jl can find r and doesnt try to install its own via conda
5. `export LD_LIBRARY_PATH=/opt/apps/anaconda/2020.07/lib:$LD_LIBRARY_PATH` # prevents cryptic errors down the line where R or julia use stuff that isnt available in the normal hpc library path
5. `export LD_LIBRARY_PATH=/opt/apps/R/4.0.4/lib64/R/lib:$LD_LIBRARY_PATH` # rcall sometimes cant find this stuff
6. `export JULIA_NO_VERIFY_HOSTS="**" && export JULIA_SSL_NO_VERIFY_HOSTS="**" && export JULIA_SSH_NO_VERIFY_HOSTS="**"` # stop julia from ever complaining about connecting to servers or downloading stuff
7. `module load zlib` # on hpc this is neccessary otherwise julia stuff that uses compression will fail
8. from within julia
    1. `] add RCall`
    2. `using Pkg`
    3. `Pkg.build("RCall")`
    4. `add git@github.com:mkarikom/DistributedTopicModels.jl.git` # install distributedtopicmodels.jl
    5. `using DistributedTopicModels` # make sure precompilation succeeds

### how to install JuliaCall (install RCall first)
1. `export JULIA_PROJECT=/dfs5/bio/mkarikom/code/DURIAN` # wherever the Project.toml with DistributedTopicModels.jl is installed
2. `export JULIA_HOME=/opt/apps/julia/1.6.0/bin` # make sure all workers can access the project enviroment
3. from within R
    1. `devtools::install_github("Non-Contradiction/JuliaCall")`

### notes on dependencies
1. The R package DrImpute calculates a distance metric based on `stats::cor` which may be unexpectedly overridden.  To prevent this, please run `devtools::install(git@github.com:mkarikom/DrImpute.git)` to install an updated fork of the package.
### To run the slurm benchmarks:
1. Please ensure that the tested versions of all julia, python, and R dependencies are installed as specified by `Mannifest.toml`,`requirements.txt`,`session_info.txt`, respectively
    * On some clusters, `~` is not guarranteed to be accessible from compute nodes during a run.  Therefore, `slurm/durian_pseudobulk/run_sim_all.sh` explicitly sets all library paths for the above dependencies under `JULIA_PROJECT`,`PYTHONPATH`,`R_LIBS_USER` and users should ensure that these are accessible.
2. Run `slurm/durian_pseudobulk/run_sim_all.sh` on the cluster, ensuring that the absolute paths for variables like `PROJECTDIR` correspond to the user.
